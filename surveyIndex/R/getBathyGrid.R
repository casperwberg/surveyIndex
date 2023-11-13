##' Get bathymetric prediction grid corresponding to the area for a DATRASraw object using the 'marmap' package  
##'
##' @title Get bathymetric prediction grid corresponding to the area for a DATRASraw object using the 'marmap' package  
##' @param d DATRASraw object
##' @param minDepth Minimum depth to include
##' @param maxDepth Maximum depth to include
##' @param resolution grid resolution (see marmap::getNOAA.bathy)
##' @param maxDist Do not include grid points farther than maxDist from nearest observation.
##' @param keep Save grid on disk for fast loading next time?
##' @param shapefile extra shapefile information to add (optional)
##' @param select columns to extract from shapefile
##' @return data.frame with depths and geographical coordinates
##' @export
getBathyGrid<-function(d,minDepth=10, maxDepth=Inf, resolution=2,maxDist=Inf, keep=TRUE, shapefile=NULL, select=NULL){
    if (!requireNamespace("marmap", quietly = TRUE)) {
        stop("Package \"marmap\" needed for this function to work. Please install it.",
      call. = FALSE)
    }
    bathy <- marmap::getNOAA.bathy(lon1=min(d$lon),lon2=max(d$lon),lat1=min(d$lat),lat2=max(d$lat),resolution=resolution, keep=keep)
    xyzbathy <- as.data.frame( marmap::as.xyz(bathy) )
    colnames(xyzbathy) <- c("lon","lat","Depth")
    xyzbathy$Depth <- -xyzbathy$Depth
    xyzbathy <- subset(xyzbathy, Depth>minDepth & Depth<maxDepth)
    if(is.finite(maxDist)){
        if (!requireNamespace("RANN", quietly = TRUE)) {
            stop("Package \"RANN\" needed for this function to work when 'maxDist' is used. Please install it.",
             call. = FALSE)
        }
        nearest <- RANN::nn2(d[[2]][,c("lon","lat")],xyzbathy[,1:2],k=1)
        xyzbathy <- xyzbathy[ nearest$nn.dists < maxDist ,]
    }
    xyzbathy <- as.data.frame(xyzbathy)
    if (is.character(shapefile)) {
        if (file.exists(shapefile)) {
            shape <- maptools::readShapeSpatial(shapefile)
        } else {
            stop("shapefile not found")
        }
        tmp <- xyzbathy
        sp::coordinates(tmp) <- ~lon + lat
        xtra <- sp::over(tmp, shape)
        if (!is.null(select)) 
            xtra <- xtra[select]
        xyzbathy <- cbind(xyzbathy, xtra)
    }
    
    xyzbathy
}
