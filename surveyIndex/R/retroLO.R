##' Make retrospective analysis for a "surveyIdx" model.
##'
##' @title Make retrospective analysis for a "surveyIdx" model.
##' @param model object of class "surveyIdx" as created by "getSurveyIdx"
##' @param d DATRASraw dataset
##' @param grid surveyIndexGrid object (see getGrid) defining the grid.
##' @param npeels number of years to successively peel of the data set
##' @return SIlist (list of surveyIndex objects)
##' @param ... Optional extra arguments to "gam"
##' @export
retro.surveyIdx<-function(model, d, grid,npeels=5,predD=NULL,...){
    if(is.null(predD)){
        predD = subset(d, haul.id %in% grid[[3]])
        predD = predD[[2]]
    }
    ages = as.numeric(colnames(model$idx))
    dataAges = model$dataAges
    famVec = model$family
    cutOff = model$cutOff

    ## Potential problem with ref gear (most common gear might change after peeling....)
    ## so add it to predfix list     
    predfix = model$predfix
    predfix$Gear = model$refGear
    
    ## How to handle partial data years?
    ## A: find last data year,quarter always remove one year from that
    lastY = max(model$yearNum)
    lastQ = max(as.numeric(as.character(d$Quarter[ d$Year == lastY ] )))
    yearRange = min(model$yearNum):max(model$yearNum)
    
    ## Recreate formulae needed
    mp <- mz <- character(length(ages))
    for(aa in 1:length(ages)){
        mp[aa] = as.character(model$pModels[[aa]]$formula)[3]
        if(length(model$zModels)>0){
            mz[aa] = as.character(model$zModels[[aa]]$formula)[3]
        } else { mz = NULL }
    }
    res = list()
    for(i in 1:npeels){
    
        curd = subset(d, Year %in% as.character(head(yearRange,length(yearRange)-(i+1))) |
                         ( Year %in% as.character(head(yearRange,length(yearRange)-i)) & 
                           Quarter %in% as.character(1:lastQ) )
                      )

        cat("Peel ",i, ": re-fitting using years ",levels(curd$Year),"\n")
        res[[i]] = getSurveyIdx(curd,ages,myids=NULL,predD=predD,cutOff=cutOff,fam=famVec,method=model$pModels[[1]]$method,knotsP=model$knotsP,knotsZ=model$knotsZ,predfix=predfix,nBoot=0,modelP=mp,modelZ=mz,...) 
    }

    class(res)<-"SIlist"
    res
    
}
##' Make leave-one-out analysis
##'
##' @title Make leave-one-out analysis
##' @param model object of class "surveyIdx" as created by "getSurveyIdx"
##' @param d DATRASraw dataset
##' @param grid surveyIndexGrid object (see getGrid) defining the grid.
##' @param fac a factor in d to leave out one at a time, e.g. d$Survey
##' @param ... Optional extra arguments to "gam"
##' @return SIlist (list of surveyIndex objects)
##' @export
leaveout.surveyIdx<-function(model,d,grid,fac,predD=NULL,...){
    stopifnot(is.factor(fac))
    
    if(is.null(predD)){
        predD = subset(d, haul.id %in% grid[[3]])
        predD = predD[[2]]
    }
    ages = as.numeric(colnames(model$idx))
    dataAges = model$dataAges
    famVec = model$family
    cutOff = model$cutOff

    predfix = model$predfix
    predfix$Gear = model$refGear

    ## Recreate formulae needed
    mp <- mz <- character(length(ages))
    for(aa in 1:length(ages)){
        mp[aa] = as.character(model$pModels[[aa]]$formula)[3]
        if(length(model$zModels)>0){
            mz[aa] = as.character(model$zModels[[aa]]$formula)[3]
        } else { mz = NULL }
    }
    
    res = list()
    
    for(facc in levels(fac)){
    
        curd = d[ which(fac!=facc) ]

        cat("Re-fitting without",facc,"\n")
        res[[facc]] = getSurveyIdx(curd,ages,myids=NULL,predD=predD,cutOff=cutOff,fam=famVec,method=model$pModels[[1]]$method,knotsP=model$knotsP,knotsZ=model$knotsZ,predfix=predfix,nBoot=0,modelP=mp,modelZ=mz,...) 
    }

    class(res)<-"SIlist"
    res
}
##' Plot survey index list (e.g. retrospective analysis)
##'
##' @title Plot survey index list (e.g. retrospective analysis)
##' @param x (named) list of "surveyIdx" objects for example from "retro.surveyIdx" or "leaveout.surveyIdx"
##' @param base Either index of x that should considered the "base run" (integer), OR object of class "surveyIdx". Confidence bounds will be shown for this model only.
##' @param rescale Should indices be rescaled to have mean 1 (over the set of intersecting years)? Default: FALSE
##' @param lwd line width argument to plot
##' @param main if not NULL override main plotting default title of "Age group a"  
##' @param allCI show 95\% confidence lines for all indices? Default FALSE.
##' @return nothing
##' @export
plot.SIlist<-function(x, base=1, rescale=FALSE,lwd=1.5,main=NULL,allCI=FALSE){
    if(class(base)=="surveyIdx"){
        x = c( list(base), x)
        base = 1
    }
    stopifnot(is.numeric(base))
    nx = length(x)
    mainwasnull = is.null(main)
    n = ncol(x[[base]]$idx)
    if(n>1){
        op <- par(mfrow=n2mfrow(n))
        on.exit(par(op))
    }

    cols = rainbow(nx) 
    if(nx==2) cols = 2:3 
    cols[base] = "black"
    allyears = lapply(x, function(x) rownames(x$idx))
    rsidx = 1:nrow(x[[nx]]$idx)
    if(rescale){
        commonyears = allyears[[1]]
        if(nx>1){
            for(i in 2:nx){
                commonyears = intersect(commonyears,allyears[[i]])
            }
            if(length(commonyears)==0) stop("rescaling not possible because the set of common years is empty")
        }
    }
    
    ss <- ssbase <- 1
        
    for(aa in 1:n){
        yl = range( c(x[[base]]$lo[,aa],x[[base]]$up[,aa],x[[base]]$idx[,aa]) )
        if(rescale){
            rsidx = which(rownames(x[[base]]$idx) %in% commonyears )
            ssbase =  mean( x[[base]]$idx[rsidx,aa], na.rm=TRUE)
            yl = yl/ssbase
        }
        if(mainwasnull) main <- paste("Age group", colnames(x[[base]]$idx)[aa])
        y = as.numeric(rownames(x[[base]]$idx))
        plot(y,x[[base]]$idx[,aa]/ssbase,type="b",ylim=yl,main=main,xlab="Year",ylab="Index")
    
        polygon(c(y, rev(y)), c(x[[base]]$lo[,aa], rev(x[[base]]$up[,aa]))/ssbase, col = "lightgrey", border = NA)
        
        for(i in 1:length(x)){
            y = as.numeric(rownames(x[[i]]$idx))
            if(rescale){
                rsidx = which(rownames(x[[i]]$idx) %in% commonyears )
                ss = mean( x[[i]]$idx[rsidx,aa], na.rm=TRUE)
            }
            lines(y,x[[i]]$idx[,aa]/ss,col=cols[i],type="b", lwd=lwd)

            if(allCI && i!=base){
                lines(y,x[[i]]$lo[,aa]/ss,col=cols[i],lwd=lwd*0.6,lty=2)
                lines(y,x[[i]]$up[,aa]/ss,col=cols[i],lwd=lwd*0.6,lty=2)
            }   
            
        }
        y = as.numeric(rownames(x[[base]]$idx))
        lines(y,x[[base]]$idx[,aa]/ssbase,type="b",lwd=lwd)
        
    }
    if(!is.null(names(x))){
        legend("topleft",legend=names(x),col=cols,lty=1,lwd=lwd,pch=1)
    }
}
##' Mohn's rho for retrospective analysis
##'
##' @title Mohn's rho for retrospective analysis
##' @param x SIlist as returned by the 'retro.surveyIdx' function
##' @param base Object of class 'surveyIdx' (full run with all years)
##' @return vector with mohn's rho for each age group
##' @export
mohn.surveyIdx<-function(x,base){
    mohns = rep(NA,ncol(base$idx))
    for(aa in 1:ncol(base$idx)){
        bias <- sapply(x,function(xx){
            y <- rownames(xx$idx)[nrow(xx$idx)]
            (xx$idx[rownames(xx$idx) == y, aa] - base$idx[rownames(base$idx) == y, aa])/base$idx[rownames(base$idx) == y,aa ]
        })
        mohns[aa] = mean(bias)
    }
  mohns   
}
