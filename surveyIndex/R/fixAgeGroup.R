##' Helper function to "borrow" missing age groups from other years
##'
##' In years where there are less than 'n' individuals of age 'age',
##' add fake individuals of that age such that there are 'n'.
##' The length of the individuals are set to the mean (or whatever 'fun' specifies)
##' of all other individuals of the same age.
##' For the minimum and maximum age groups fun it is reasonable to replace 'mean' with 'min' and 'max' respectively.
##' Note, that you might need to call 'addSpectrum' on the object again.
##' @title Helper function to "borrow" missing age groups from other years
##' @param x DATRASraw object
##' @param age age to impute
##' @param n at least this many individuals in each year
##' @param fun A function such as 'mean','median','min', or 'max'.
##' @return a DATRASraw object
##' @export
fixAgeGroup <-
function(x,age=0,n=3,fun="mean"){
  cm.breaks<-attr(x,"cm.breaks")  
  f <- match.fun(fun)
  d=split(x,x$Year)
  subsLength=f(x[[1]]$LngtCm[x[[1]]$Age==age],na.rm=TRUE)
  for(y in 1:length(d)){
        nobs = sum(d[[y]][[1]]$Age==age,na.rm=TRUE)
        if(nobs<n) {
            sel=sample(1:nrow(d[[y]][[1]]),n-nobs)
            d[[y]][[1]]=rbind(d[[y]][[1]][sel,],d[[y]][[1]]);
            d[[y]][[1]][1:(n-nobs),"Age"]=age;
            d[[y]][[1]][1:(n-nobs),"LngtCm"]=subsLength;
            d[[y]][[1]][1:(n-nobs),"NoAtALK"]=1;
        }
      }
  dd <- do.call("c",d)
  if(!is.null(cm.breaks)) dd<-addSpectrum(dd,cm.breaks=cm.breaks)
  dd
}
