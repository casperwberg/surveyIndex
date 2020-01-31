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
##' @param x object of class "SIlist" as created by e.g. "retro.surveyIdx" or "leaveout.surveyIdx"
##' @param base object of class "surveyIdx" (base model)
##' @param rescale Should indices be rescaled to have mean 1 (over the set of intersecting years)? Default: FALSE
##' @param lwd line width argument to plot
##' @return nothing
##' @export
plot.SIlist<-function(x, base=x[[1]], rescale=FALSE,lwd=1.5,main=NULL, basename = "base"){
    mainwasnull = is.null(main)
    n = ncol(base$idx)
    op <- par(mfrow=n2mfrow(n))
    on.exit(par(op))
    nx = length(x)
    cols = rainbow(nx)
    allyears = lapply(x, function(x) rownames(x$idx))
    rsidx = 1:nrow(x[[nx]]$idx)
    if(rescale && nx>1){
        commonyears = allyears[[1]]
        for(i in 2:nx){
            commonyears = intersect(commonyears,allyears[[i]])
        }
        if(length(commonyears)==0) stop("rescaling not possible because the set of common years is empty")
    }
    
    ss = 1
        
    for(aa in 1:n){
        yl = range( c(base$lo[,aa],base$up[,aa],base$idx[,aa]) )
        if(rescale){
            rsidx = which(rownames(base$idx) %in% commonyears )
            ss =  mean( base$idx[rsidx,aa], na.rm=TRUE)
            yl = yl/ss
        }
        if(mainwasnull) main <- paste("age group", colnames(base)[aa])
        y = as.numeric(rownames(base$idx))
        plot(y,base$idx[,aa]/ss,type="b",ylim=yl,main=main,xlab="Year",ylab="Index")
    
        polygon(c(y, rev(y)), c(base$lo[,aa], rev(base$up[,aa]))/ss, col = "lightgrey", border = NA)
        lines(y,base$idx[,aa]/ss,type="b",lwd=lwd)

        for(i in 1:length(x)){
            y = as.numeric(rownames(x[[i]]$idx))
            if(rescale){
                rsidx = which(rownames(x[[i]]$idx) %in% commonyears )
                ss = mean( x[[i]]$idx[rsidx,aa], na.rm=TRUE)
            }
            lines(y,x[[i]]$idx[,aa]/ss,col=cols[i],type="b", lwd=lwd)
        }
        
    }
    if(!is.null(names(x))){
        legend("topleft",legend=c(basename,names(x)),col=c("black",cols),lty=1,lwd=lwd,pch=1)
    }
}
