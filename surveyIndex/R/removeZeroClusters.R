##' Remove factor levels with only/many zero observations, which can cause estimation problems
##'
##' @title Remove factor levels with only zero observations, which can cause estimation problems
##' @param d DATRASraw object
##' @param cutat numeric in the interval [0,1[. If greater than 0, then factor levels where mu/max(mu) < cutat are removed. Note, that this will remove some positive observations.   
##' @param response Character naming the response column. If NULL, then total numbers are used as response.
##' @param factors vector naming the factors to consider.
##' @param verbose print out extra info? Default: TRUE
##' @return DATRASraw object
##' @export
removeZeroClusters<-function(d,cutat=0,response=NULL,factors=c("Gear","Ship","Year","StatRec"),verbose=TRUE){

    N = nrow(d[[2]])
    if(is.null(response)){
        if(is.null(attr(d, "cm.breaks"))) suppressWarnings(d <- addSpectrum(d))
        d$Ntot = rowSums(d$N)
        response = "Ntot"
    }

    ## Remove empty/all zero levels
    for(fac in factors){
        factab <- xtabs(d[[2]][,response] ~ d[[2]][,fac])
        
        if(any(factab==0)){
            if(verbose) cat("Removing following empty/all zero levels of ", fac,":",names(factab[factab==0]),"\n")
            badids <- d$haul.id[ d[[2]][,fac] %in% names(factab[factab==0]) ]
            d <- subset(d, !haul.id %in% badids) 
        }

        ## remove also levels where mean catch is less than cutat*max(mean) 
        factab <- xtabs(d[[2]][,response] ~ d[[2]][,fac])
        ntab <- xtabs( ~ d[[2]][,fac])
        mutab <- factab / ntab
        maxmu <- max(mutab)

        nposremoved <- 0
        if(cutat>0 & any(mutab/maxmu < cutat)){
            if(verbose) cat("Removing following levels of ", fac,":",names(factab[mutab/maxmu < cutat]),"\n")
            badids <- d$haul.id[ d[[2]][,fac] %in% names(factab[mutab/maxmu < cutat]) ]
            bad <- subset(d, haul.id %in% badids)
            d <- subset(d, !haul.id %in% badids)
            nposremoved <- sum(bad[[2]][,response]>0)
        }
        
        ## re-factor (drop empty factor levels)
        d[[2]][,fac] = factor(d[[2]][,fac],levels=unique(d[[2]][,fac]))
    }

    Nnew = nrow(d[[2]])
    if(verbose) {
        cat("Number of hauls before: " , N, " , after: ",Nnew, " (",round((N-Nnew)/N*100,1),"% reduction)\n")
        cat("Number of positive hauls removed: ",nposremoved,"\n")
    }
    
    return(d)
    
}
