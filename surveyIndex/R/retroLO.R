##' Make retrospective analysis for a "surveyIdx" model.
##'
##' @title Make retrospective analysis for a "surveyIdx" model.
##' @param model object of class "surveyIdx" as created by "getSurveyIdx"
##' @param d DATRASraw dataset
##' @param grid surveyIndexGrid object (see getGrid) defining the grid.
##' @param npeels number of years to successively peel of the data set
##' @return SIlist (list of surveyIndex objects)
##' @export
retro.surveyIdx<-function(model, d, grid, npeels=5){
    predD = subset(d, haul.id %in% grid[[3]])
    predD = predD[[2]]
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
    
    ## Problem if grid includes points from "npeels" recent years...
    ## Solution: use "predD" for the grid
    res = list()
    for(i in 1:npeels){
    
        curd = subset(d, Year %in% as.character(head(yearRange,length(yearRange)-(i+1))) |
                         ( Year %in% as.character(head(yearRange,length(yearRange)-i)) & 
                           Quarter %in% as.character(1:lastQ) )
                      )

        cat("Peel ",i, ": re-fitting using years ",levels(curd$Year),"\n")
        res[[i]] = getSurveyIdx(curd,ages,myids=NULL,predD=predD,cutOff=cutOff,fam=famVec,method=model$pModels[[1]]$method,knotsP=model$knotsP,knotsZ=model$knotsZ,predfix=predfix,nBoot=0) 
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
##' @return SIlist (list of surveyIndex objects)
##' @export
leaveout.surveyIdx<-function(model,d,grid,fac){
    predD = subset(d, haul.id %in% grid[[3]])
    predD = predD[[2]]
    ages = as.numeric(colnames(model$idx))
    dataAges = model$dataAges
    famVec = model$family
    cutOff = model$cutOff

    predfix = model$predfix
    predfix$Gear = model$refGear
    res = list()
    
    for(facc in levels(fac)){
    
        curd = d[ which(fac!=facc) ]

        cat("Re-fitting without",facc,"\n")
        res[[facc]] = getSurveyIdx(curd,ages,myids=NULL,predD=predD,cutOff=cutOff,fam=famVec,method=model$pModels[[1]]$method,knotsP=model$knotsP,knotsZ=model$knotsZ,predfix=predfix,nBoot=0) 
    }

    class(res)<-"SIlist"
    res
}

##' Plot survey index list (e.g. retrospective analysis)
##'
##' @title Plot survey index list (e.g. retrospective analysis)
##' @param x object of class "SIlist" as created by e.g. "retro.surveyIdx" or "leaveout.surveyIdx"
##' @param base 
##' @return nothing
##' @export
plot.SIlist<-function(x, base=x[[1]],lwd=1.5){

    n = ncol(base$idx)
    op <- par(mfrow=n2mfrow(n))
    on.exit(par(op))
    nx = length(x)
    cols = rainbow(nx)
    
    for(aa in 1:n){
        y = as.numeric(rownames(base$idx))
        yl = range( c(base$lo[,aa],base$up[,aa],base$idx[,aa]) )
        plot(y,base$idx[,aa],type="b",ylim=yl,main=paste0("age group ",aa),xlab="Year",ylab="Index")
    
        polygon(c(y, rev(y)), c(base$lo[,aa], rev(base$up[,aa])), col = "lightgrey", border = NA)
        lines(y,base$idx[,aa],type="b",lwd=lwd)

        for(i in 1:length(x)){
            y = as.numeric(rownames(x[[i]]$idx))
            lines(y,x[[i]]$idx[,aa],col=cols[i],type="b", lwd=lwd)
        }
        
    }
    if(!is.null(names(x))){
        legend("topleft",legend=names(x),col=cols,lty=1,lwd=lwd,pch=1)
    }
}
