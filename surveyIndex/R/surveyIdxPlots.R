##' Visualize results from a survey index model fitted with getSurveyIdx().
##'
##' @title Visualize results from a survey index model fitted with getSurveyIdx().
##' @param x Survey index as produced by getSurveyIndex()
##' @param dat DATRASraw object
##' @param alt.idx optional matrix with alternative index
##' @param myids vector of haul ids that constitute the grid
##' @param cols which age columns to consider?
##' @param select character vector of chosen plots. Either one of "index","map","residuals", or "fitVsRes" or a number. Numbers refer to smooths in the order they appear in the formula.
##' @param par 'par' settings for plotting (a named list).
##' @param colors colors for spatial effect.
##' @param map.cex size of grid points on maps
##' @param plotByAge boolean (default=TRUE). If true, par(par) is called for each age group.
##' @param legend boolean (default=TRUE). add legends to plot?
##' @param predD DATRASraw object with grid (optional). Overrides 'myids' if supplied.
##' @param year numeric (default=NULL). If 'select' equals 'map' a specific year can be chosen (only meaningful for time-varying spatial effects). 
##' @param ... Additional parameters for plot()
##' @return nothing
##' @export
##' @import maps mapdata tweedie
surveyIdxPlots <-
function(x,dat,alt.idx=NULL,myids,cols=1:length(x$pModels),select=c("index","map","residuals","fitVsRes"),par=list(mfrow=c(3,3)),colors=rev(gray.colors(5)),map.cex=1,plotByAge=TRUE,legend=TRUE,predD=NULL,year=NULL,...){

  if(!plotByAge & !is.null(par)) par(par)
  for(a in cols){
    if(plotByAge & !is.null(par) ) par(par)

    if(any(select=="index")){
      ys=range(as.numeric(levels(dat$Year)));
      ys=ys[1]:ys[2]
      yl=range( c(x$idx[,a]/mean(x$idx[,a],0) ))*1.10
      if(!is.null(alt.idx) && a<=ncol(alt.idx)) {
          yl=range( c(alt.idx[,a]/mean(alt.idx[,a]),x$idx[,a]/mean(x$idx[,a],0) ) )*1.10
          plot(ys,alt.idx[,a]/mean(alt.idx[,a],na.rm=TRUE),ylim=yl,col=2,ylab="Index",xlab="Year",main=paste("Age group",a)) } else {
        plot(ys,rep(NA,length(ys)),ylim=yl,col=2,ylab="Index",xlab="Year",main=paste("Age group",a))
      }

      idx=x$idx
      lo=x$lo
      up=x$up
      idx[x$idx<=0]=NA
      lo[x$idx<=0]=NA
      up[x$idx<=0]=NA 
      
      lines(ys,idx[,a]/mean(idx[,a],na.rm=TRUE),lwd=2)
      lines(ys,lo[,a]/mean(idx[,a],na.rm=TRUE),lwd=2,lty=2);
      lines(ys,up[,a]/mean(idx[,a],na.rm=TRUE),lwd=2,lty=2);
      
      if(legend) legend("topright",pch=c(1,NA),lty=c(NA,1),col=c(2,1),legend=c("alt.idx","GAM"));
    }

    if(any(select=="map")){
        
        xlims=range(dat$lon,na.rm=TRUE)
        ylims=range(dat$lat,na.rm=TRUE)
      if(is.null(predD)){ tmp=subset(dat,haul.id %in% myids) } else {tmp=predD;}
      
      if(is.null(year)) { concT = concTransform(log(x$gPreds[[a]])) }  else {
                           y = which( as.numeric(as.character(names(x$gPreds2[[a]])))==year)
                           if(length(y)==0) stop(paste("Year",year, "age group",a,"not found."))
                           concT = concTransform(log(x$gPreds2[[a]][[y]]))
                     }
      if(length(colors)>1) zFac=cut( concT,0:length(colors)/length(colors)) else zFac=1;
      if(length(map.cex)>1) sFac=cut(log(x$gPreds[[a]]),length(map.cex)) else sFac=1;
      myCols=colors;
      plot(tmp$lon,y=tmp$lat,col=1,pch=1,cex=map.cex[sFac],xlim=xlims,ylim=ylims,xlab="Longitude",ylab="Latitude",main=paste("Age group",a,year),...)
      points(tmp$lon,y=tmp$lat,col=myCols[zFac],pch=16,cex=map.cex[sFac]-0.05)
      maps::map('worldHires',xlim=xlims,ylim=ylims,fill=TRUE,plot=TRUE,add=TRUE,col=grey(0.5))
        
      if(legend) legend("topright",legend=levels(zFac),pch=16,col=colors,bg="white")
    }
    for(k in 1:length(select)){
        ss = suppressWarnings(as.numeric(select[k]));
        if(!is.na(ss)){
            plot.gam(x$pModels[[a]],select=ss,main=paste("age gr",a),...);
        }
    }

    if(any(select=="residuals") || any(select=="fitVsRes") || any(select=="resVsYear") || any(select=="resVsShip") || any(select=="spatialResiduals")){
        if(pmatch("Tweedie",x$pModels[[a]]$family$family,nomatch=-1)==1){
            resi <- qres.tweedie(x$pModels[[a]]);
        } else {
            resi <- residuals(x$pModels[[a]]);
        }
    }
    
    if(any(select=="residuals")){
        hist(resi,nclass=30,main=paste("Residuals (pos only) age gr",a),xlab="Residuals")
    }
    if(any(select=="fitVsRes")){
      plot(fitted(x$pModels[[a]]),resi,xlab="Fitted",ylab="Residuals",main=paste("age gr",a),...)
    }

    if(any(select=="resVsYear")){
       plot(x$pData[[a]]$Year,resi,main=paste("age gr",a),xlab="Year",ylab="Residuals",...)
    }

    if(any(select=="resVsShip")){
       plot(x$pData[[a]]$Ship,resi,main=paste("age gr",a),xlab="Year",ylab="Residuals",...)
    }

    if(any(select=="spatialResiduals")){
  
        scale <- 3
        if(is.null(year)) stop("a year must be supplied")
        sel <- which(x$pData[[a]]$Year == as.character(year))
        plot(x$pData[[a]]$lon, x$pData[[a]]$lat, type = "n", xlab = "Longitude", ylab = "Latitude",main=paste("Age group",a,year))
        maps::map("worldHires", fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
        
        positive = resi>0
        points(x$pData[[a]]$lon[sel][positive], x$pData[[a]]$lat[sel][positive], pch = 1, cex = scale * sqrt(resi[positive]),col="blue")
        points(x$pData[[a]]$lon[sel][!positive], x$pData[[a]]$lat[sel][!positive], pch = 1, cex = scale * sqrt(-resi[!positive]),col="red")
        
    }
    
  }

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
