##' Likelihood ratio test for comparing two survey indices.
##' 
##' @title Likelihood ratio test for comparing two survey indices.
##' @param m1 
##' @param m2 
##' @return A p-value.
##' @export anova.SI
anova.SI <-
function(m1,m2){
  ll1=m1$logLik
  ll2=m2$logLik
  edfs1=m1$edfs
  edfs2=m2$edfs
  
  if (edfs2 > edfs1) {
    1 - pchisq(2 * (ll2 - ll1), edfs2 - edfs1)
  }
  else {
    1 - pchisq(2 * (ll1 - ll2), edfs1 - edfs2)
  }
}

##' Akaike Information Criterion (or BIC) for survey index models
##'
##' @title Akaike Information Criterion (or BIC) for survey index models
##' @param x survey index as return from getSurveyIdx
##' @param BIC if TRUE compute BIC instead of AIC
##' @return numeric value
##' @export AIC.surveyIdx
AIC.surveyIdx<-function(x, BIC=FALSE){
    if(!BIC) return(2*x$edfs-2*x$logLik)
    if(pmatch("Tweedie",x$pModels[[1]]$family$family,nomatch=-1)==1 ||
       pmatch("negbin",x$pModels[[1]]$family$family,nomatch=-1)==1 ){
        log(length(x$pModels[[1]]$y)*(length(x$pModels)))*x$edfs-2*x$logLik
    } else {
             log(length(x$zModels[[1]]$y)*(length(x$zModels)))*x$edfs-2*x$logLik
    }
}
