##' Create a grid of haul positions from a DATRASraw object.
##'
##' @title Create a grid of haul positions from a DATRASraw object.
##' @param dd DATRASraw object
##' @param nLon number of grid cells in the longitude direction.
##' @return A surveyIndexGrid (a list of coordinates and haul.ids)
##' @export
getGrid <-
function(dd,nLon=20){
  mlon=mean(dd$lon)
  mlat=mean(dd$lat)
  kmPerDegLon=gcd.hf(deg2rad(mlon),deg2rad(mlat),deg2rad(mlon+1),deg2rad(mlat))
  kmPerDegLat=gcd.hf(deg2rad(mlon),deg2rad(mlat),deg2rad(mlon),deg2rad(mlat+1))

  getBps <- function(labs){
        cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
              upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  }

  lonf = cut(dd$lon, nLon,dig.lab=8)
  lonbps = getBps(levels(lonf))
  gridsize = diff(lonbps[1,])*kmPerDegLon
  nLat = round(diff(range(dd$lat))*kmPerDegLat/gridsize)
  latf = cut(dd$lat, nLat, ,dig.lab=8)
  dd$StatRec2 = as.factor(paste(lonf,latf))

  latbps=getBps(levels(latf))
  cat("Approximate grid size: ", mean( c(diff(lonbps[1,])*kmPerDegLon,diff(latbps[1,])*kmPerDegLat))," km\n");
  
  uRecs=unique(as.character(dd$StatRec2))
  N=length(uRecs);
  mylon=numeric(N);
  mylat=numeric(N);
  myids=numeric(N);
  k=0;
  for(rec in uRecs)
    {
      k=k+1;
      tmp=subset(dd[[2]],StatRec2==rec);
      mlon=mean(tmp$lon);
      mlat=mean(tmp$lat);
      dist=sqrt( (mlon-tmp$lon)^2+(mlat-tmp$lat)^2);
      sel=which.min(dist);
      mylon[k]=tmp$lon[sel];
      mylat[k]=tmp$lat[sel];
      myids[k]=as.character(tmp$haul.id[sel]);
    }
  ret <- list(mylon,mylat,myids,lonbps,latbps)
  class(ret) <- "surveyIndexGrid"
  return(ret);
}

##' Plot a surveyIndexGrid 
##'
##' @title Plot a surveyIndexGrid 
##' @param grid a surveyIndexGrid (as created by the "getGrid" function)
##' @return nothing
##' @export
plot.surveyIndexGrid<-function(grid, pch=1,gridCol="lightgrey"){
    plot(grid[[1]],grid[[2]],xlab="Longitude",ylab="Latitude",pch=pch)
    lonbps = grid[[4]]
    latbps = grid[[5]]
    for(i in 1:nrow(lonbps)) abline(v=lonbps[i,1],col=gridCol)
    abline(v=tail(lonbps,1)[2],col=gridCol)
    for(i in 1:nrow(latbps)) abline(h=latbps[i,1],col=gridCol)
    abline(h=tail(latbps,1)[2],col=gridCol)
    maps::map("worldHires", fill = TRUE, plot = TRUE,
              add = TRUE, col = grey(0.5))
}
