##' Randomized quantile residuals for class 'surveyIndex'
##'
##' @title Randomized quantile residuals for class 'surveyIndex'
##' @param x An object of type 'surveyIndex' as created by 'getSurveyIdx'
##' @param a age group
##' @return A vector of residuals, which should be iid standard normal distributed
##' @export
residuals.surveyIdx <- function(x,a=1){
    if (pmatch("Tweedie", x$pModels[[a]]$family$family, 
               nomatch = -1) == 1) {
        resi = qres.tweedie(x$pModels[[a]])
    } else if(pmatch("LogNormal", x$pModels[[a]]$family$family, 
               nomatch = -1) == 1 ){## Delta-Lognormal
        y = x$zModels[[a]]$y
        p0 = predict(x$zModels[[a]],type="response")
        vari = x$family$variance
        
        cdfpos = pnorm( q=y,mean=predict(x$pModels[[a]],type="link"), sd=sqrt(vari))
        u = (1-p0) + cdfpos*p0
        u[ y == 0 ] = runif(sum(y == 0), min = 0, max = u[y == 0])
        resi = qnorm(u)
            
    } else if(pmatch("Gamma", x$pModels[[a]]$family$family, 
                     nomatch = -1) == 1 ){
        requireNamespace("MASS")
        y = x$zModels[[a]]$y
        p0 = predict(x$zModels[[a]],type="response")
        means = predict(x$pModels[[a]],type="response")
        shape = MASS::gamma.shape(x$pModels[[a]])$alpha
        cdfpos = pgamma( q=y,shape=shape,rate=shape/means) 
        u = (1-p0) + cdfpos*p0
        u[ y == 0 ] = runif(sum(y == 0), min = 0, max = u[y == 0])
        resi = qnorm(u)
        
    } else if(pmatch("Negative Binomial", x$pModels[[a]]$family$family, 
                     nomatch = -1) == 1 ){
        y = x$pModels[[a]]$y
        size = x$pModels[[a]]$family$getTheta(TRUE)
        mu = fitted(x$pModels[[a]])
        p = size/(mu + size)
        a = ifelse(y > 0, pbeta(p, size, pmax(y, 1)), 0)
        b = pbeta(p, size, y + 1)
        u = runif(n = length(y), min = a, max = b)
        resi = qnorm(u)
    }
    
    resi
}

##' Randomized quantile residuals for Tweedie models
##'
##' @title Randomized quantile residuals for Tweedie models
##' @param gam.obj A gam object (mgcv package)
##' @return A vector of residuals, which should be iid standard normal distributed
##' @export
##' @import tweedie
qres.tweedie<-function (gam.obj) 
{
    requireNamespace("tweedie")
    mu <- fitted(gam.obj)
    y <- gam.obj$y
    df <- gam.obj$df.residual
    w <- gam.obj$prior.weights
    if (is.null(w)) 
        w <- 1
    p <- gam.obj$family$getTheta(TRUE)
    dispersion <- gam.obj$scale
    u <- tweedie::ptweedie(q = y, power = p, mu = fitted(gam.obj), 
        phi = dispersion/w)
    if (p > 1 && p < 2) 
        u[y == 0] <- runif(sum(y == 0), min = 0, max = u[y == 
            0])
    qnorm(u)
}
