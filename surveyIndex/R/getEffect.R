##' Calculate confidence intervals for a named parameter in a survey index model.
##'
##' @title Calculate confidence intervals for a named parameter in a survey index model.
##' @param x survey index
##' @param dat DATRASraw object
##' @param parName name of the parameter, e.g. "Gear"
##' @param cutOff see getSurveyIndex()
##' @param nboot see getSurveyIndex()
##' @param pOnly only calculate for positive part of model, defaults to FALSE.
##' @return list of estimates + ci bounds for each age group.
##' @importFrom MASS mvrnorm
##' @export
getEffect <-
function(x,dat,parName="Gear",cutOff,nboot=1000,pOnly=FALSE){
    noAges=length(x$pModels);

    res=list()
    for(a in 1:noAges){
        cat("Age ",a,"\n");
        shipSelP = grep(parName,names(coef(x$pModels[[a]])))
        shipSelZ = grep(parName,names(coef(x$zModels[[a]])))

        dd=subset(dat,Nage[,a]>cutOff)
        zNam=levels(dat$Ship)
        pNam=levels(dd$Ship)
        dif=setdiff(zNam,pNam);
        remo = which(zNam %in% dif)
        if(length(remo)>0) shipSelZ = shipSelZ[-remo];
        
        if(length(shipSelP)!=length(shipSelZ)) { print("unequal number of ship effects"); }

        
        brp.1=mvrnorm(n=nboot,coef(x$pModels[[a]]),x$pModels[[a]]$Vp);
        brp.0=mvrnorm(n=nboot,coef(x$zModels[[a]]),x$zModels[[a]]$Vp);
        ilogit<-function(x) 1/(1+exp(-x));

        shipE = exp(brp.1[,shipSelP,drop=FALSE]);
        if(!pOnly) shipE = shipE*ilogit(brp.0[,shipSelZ,drop=FALSE])

        upres  = apply(shipE,2, quantile,probs=0.975);
        lores  = apply(shipE,2, quantile,probs=0.025);

        tmp=cbind(colMeans(shipE),upres,lores);
        
        res[[a]]=tmp;
    }
    return(res);
}
