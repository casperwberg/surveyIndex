##' Write survey index to file in standard XSA/SAM format
##'
##' @title Write survey index to file in standard XSA/SAM format
##' @param x matrix with survey indices
##' @param ages vector of ages
##' @param years vector of years
##' @param toy fraction of year the survey is conducted (between 0 and 1)
##' @param file filename to write to
##' @param nam file description header
##' @return nothing
##' @export
exportSI <-
function(x,ages,years,toy,file,nam="") {
  cat(nam,"\n",file=file)
  cat(range(as.numeric(as.character(years))),"\n",file=file,append=TRUE)
  cat("1 1 ",rep(toy,2),"\n",file=file,append=TRUE)
  cat(min(ages),max(ages),"\n",file=file,append=TRUE)
  colstokeep <- colnames(x) %in% ages
  rowstokeep <- rownames(x) %in% years
  if (all(isFALSE(colstokeep))) stop("None of the selected ages are in `x`")
  if (all(isFALSE(rowstokeep))) stop("None of the selected years are in `x`")
  write.table(round(cbind(1,x[rowstokeep,colstokeep]),4),file=file,row.names=FALSE,col.names=FALSE,append=TRUE)
}
