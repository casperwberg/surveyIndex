##' Likelihood ratio test for comparing two survey indices.
##' 
##' @title Likelihood ratio test for comparing two survey indices.
##' @param m1 
##' @param m2 
##' @return A p-value.
##' @export
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
