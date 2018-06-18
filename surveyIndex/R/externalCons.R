##' Calculate external consistencies between two survey indices.
##'
##' Proper alignment of years and ages must be ensured by the user.
##' @title Calculate external consistencies between two survey indices.
##' @param tt A matrix with survey indices (rows=years, cols=ages)
##' @param tt2 A matrix with survey indices (rows=years, cols=ages)
##' @param do.plot plot it?
##' @return A vector of correlations (consistencies)
##' @export
externalCons <-
function(tt,tt2,do.plot=FALSE){
  tt[tt==0]=NA
  tt2[tt2==0]=NA
  if(do.plot){ X11(); b=ceiling((ncol(tt))/2); par(mfrow=c(b,b));}
  for(a in 1:ncol(tt)){
    cat("Survey 1 Age ",a," vs Survey 2 ",a," : ",cor(log(tt[,a]),log(tt2[,a]),use="pairwise.complete.obs"),"\n")
    if(do.plot) {plot(log(tt[,a]),log(tt2[,a])); abline(0,1);}
  }
  return( sapply(1:ncol(tt),function(a) cor(log(tt[,a]),log(tt2[,a]),use="pairwise.complete.obs") ));
}
