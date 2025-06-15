
##' Survey index using the stratified mean method using ICES statistical rectangles as strata.
##'
##' @title Survey index using the stratified mean method using ICES statistical rectangles as strata.
##' @param x DATRASraw object. Must contain a matrix: x[[2]]$Nage.
##' @param ageCols which columns of the Nage matrix should be included?
##' @param doLog log-transform?
##' @return a matrix with survey indices
##' @export
getSurveyIdxStratMean <-
function(x,ageCols,doLog=FALSE){
  ysplit=split(x,x$Year);
  res=matrix(NA,nrow=length(ysplit),ncol=length(ageCols))
  for(y in 1:length(ysplit)){

    if(!doLog) {
      byRec=aggregate(ysplit[[y]][[2]]$Nage[,ageCols],by=list(ysplit[[y]][[2]]$StatRec),FUN="mean")  } else {
        byRec=aggregate(log(ysplit[[y]][[2]]$Nage[,ageCols]+1),by=list(ysplit[[y]][[2]]$StatRec),FUN="mean")
      }
    
    if(length(ageCols) == 1){
      res[y,] = mean(byRec[,-1])
    }else res[y,]=colMeans(byRec[,-1])
    }
  res
}
