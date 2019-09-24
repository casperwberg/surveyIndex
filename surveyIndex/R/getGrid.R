##' Create a grid of haul positions from a DATRASraw object.
##'
##' @title Create a grid of haul positions from a DATRASraw object.
##' @param dd DATRASraw object
##' @param nLon number of grid cells in the longitude direction.
##' @return a list of coordinates and haul.ids.
##' @export
getGrid <-
function(dd,nLon=20){
  mlon=mean(dd$lon)
  mlat=mean(dd$lat)
  kmPerDegLon=gcd.hf(deg2rad(mlon),deg2rad(mlat),deg2rad(mlon+1),deg2rad(mlat))
  kmPerDegLat=gcd.hf(deg2rad(mlon),deg2rad(mlat),deg2rad(mlon),deg2rad(mlat+1))
  nLat=round(nLon*(kmPerDegLat/kmPerDegLon))

  dd$StatRec2=as.factor(paste(cut(dd$lon,nLon),cut(dd$lat,nLat)));

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
  return(list(mylon,mylat,myids));
}
