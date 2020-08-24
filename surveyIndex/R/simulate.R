##' Simulate data from a surveyIdx model (experimental and subject to change)
##'
##' @title Simulate data from a surveyIndex model (experimental and subject to change)
##' @param model object of class 'surveyIdx' 
##' @param d A dataset (DATRASraw object) 
##' @param sampleFit Use a random sample from the gaussian approximation to distribution of the estimated parameter vector. Default: FALSE. 
##' @param condSim optional results of previous call to this function. Use this if you want to generate many datasets (much faster, since mean predictions are re-used). 
##' @return list with  1) simulated observations with noise 2) mean (no noise) 3) zero probability.
##' @export
surveyIdx.simulate<-function(model,d,sampleFit=FALSE,condSim=NULL){
    ages = as.numeric(colnames(model$idx))
    dataAges <- model$dataAges
    famVec = model$family

    out = d$Nage
    out.mu = list()
    out.mu0 = list()
    
    for(a in 1:length(ages)){


        if(sampleFit && is.null(condSim)){
            m.pos = model$pModels[[a]]
            
            Xp.1=predict(m.pos,newdata=d[[2]],type="lpmatrix");
            brp.1=MASS::mvrnorm(n=1,coef(m.pos),m.pos$Vp);
            OS.pos = matrix(0,nrow(d[[2]]),1);
            terms.pos=terms(m.pos)
            if(!is.null(m.pos$offset)){
                off.num.pos <- attr(terms.pos, "offset") 
                for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos, 
                                                                    "variables")[[i + 1]], d[[2]])
            }
            mu = Xp.1%*%brp.1+OS.pos
            out.mu[[a]] = exp(mu) ## obs, assumes log link!
            
            if(!famVec[a] %in% c("Tweedie","negbin")){
                m0 = model$zModels[[a]]
                Xp.0=predict(m0,newdata=d[[2]],type="lpmatrix");
                brp.0=MASS::mvrnorm(n=1,coef(m0),m0$Vp);
                OS0 = matrix(0,nrow(d[[2]]),1);
                terms.0=terms(m0)
                if(!is.null(m0$offset)){
                    off.num.0 <- attr(terms.0, "offset")
                    for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0, 
                                                                "variables")[[i + 1]], d[[2]])
                }
                mu0=m0$family$linkinv(Xp.0%*%brp.0+OS0);
                out.mu0[[a]] = mu0
            }
            
        }
        if(!is.null(condSim)) out.mu[[a]] = condSim[[2]][[a]]
        if(!is.null(condSim) && !famVec[a] %in% c("Tweedie","negbin")) out.mu0[[a]] = condSim[[3]][[a]]
        
        if(famVec[a]==c("Tweedie")){
            p = model$pModels[[a]]$family$getTheta(TRUE)
            phi = model$pModels[[a]]$scale
            if(!sampleFit && is.null(condSim)) {  out.mu[[a]] = predict(model$pModels[[a]],newdata=d[[2]],type="response") } 
            out[,a] = rTweedie(out.mu[[a]],p,phi)
        } else if(famVec[a]=="LogNormal"){
            if(!sampleFit && is.null(condSim)) out.mu0[[a]] = predict(model$zModels[[a]],newdata=d[[2]],type="response")
            pa = rbinom(nrow(out),1,prob=out.mu0[[a]])
            sig2 = model$pModels[[a]]$sig2
            if(!sampleFit && is.null(condSim)) out.mu[[a]] = exp(predict(model$pModels[[a]],newdata=d[[2]]))
            pos = exp( rnorm(nrow(out),log(out.mu[[a]]),sqrt(sig2)) )
            out[,a] = pa*pos
        } else if(famVec[a]=="Gamma"){
            if(!sampleFit && is.null(condSim)) out.mu0[[a]] = predict(model$zModels[[a]],newdata=d[[2]],type="response")
            pa = rbinom(nrow(out),1,prob=out.mu0[[a]])
            
            if(!sampleFit && is.null(condSim)){
                out.mu[[a]]= predict(model$pModels[[a]],newdata=d[[2]],type="response")
            }             
            shape = 1/model$pModels[[a]]$scale
            pos = rgamma(nrow(out),rate=shape/out.mu[[a]],shape=shape)
            out[,a] = pa*pos
        }

    }
    list(sim=out,mu=out.mu,prob0=out.mu0)
}
