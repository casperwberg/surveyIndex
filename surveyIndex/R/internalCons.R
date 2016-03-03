##' Calculate internal consistency of a survey index.
##'
##' @title Calculate internal consistency of a survey index.
##' @param tt A matrix with survey indices (rows=years, cols=ages)
##' @param do.plot Plot it?
##' @return a vector of consistencies
##' @export
internalCons <-
function(tt,do.plot=FALSE){
  tt[tt==0]=NA
  sel1=1:(nrow(tt)-1);
  sel2=2:nrow(tt);
  if(do.plot){ X11(); b=ceiling((ncol(tt)-1)/2); par(mfrow=c(b,b));}
  for(a in 1:(ncol(tt)-1)){
    cat("Age ",a," vs ",a+1," : ",cor(log(tt[sel1,a]),log(tt[sel2,a+1]),use="pairwise.complete.obs"),"\n")
    if(do.plot) {plot(log(tt[sel1,a]),log(tt[sel2,a+1])); abline(0,1);}
  }
  return( sapply(1:(ncol(tt)-1), function(a) cor(log(tt[sel1,a]),log(tt[sel2,a+1]),use="pairwise.complete.obs")));
}
