##' Visualize results from a survey index model fitted with getSurveyIdx().
##'
##' @title Visualize results from a survey index model fitted with getSurveyIdx().
##' @param x Survey index as produced by getSurveyIndex()
##' @param dat DATRASraw object
##' @param alt.idx optional matrix with alternative index
##' @param myids vector of haul ids that constitute the grid
##' @param cols which age columns to consider?
##' @param select character vector of chosen plots. Either one of "index","map","absolutemap","CVmap","residuals","fitVsRes",""resVsYear","resVsShip","spatialResiduals", or a number. Numbers refer to smooths in the order they appear in the formula.
##' @param par 'par' settings for plotting (a named list).
##' @param colors colors for spatial effect.
##' @param map.cex size of grid points on maps
##' @param plotByAge boolean (default=TRUE). If true, par(par) is called for each age group.
##' @param legend boolean (default=TRUE). add legends to plot?
##' @param predD DATRASraw object with grid (optional). Overrides 'myids' if supplied.
##' @param year numeric scalar or vector (default=NULL). If 'select' equals 'map' a specific year can be chosen (only meaningful for time-varying spatial effects). If select equals 'absolutemap' or 'CVmap' then year must be a vector. 
##' @param main optional main title (overrides default title)
##' @param legend.signif Number of significant digits in map legends
##' @param legend.pos Position of legend (e.g. "bottomleft") see ?legend
##' @param restoreOldPar restore old par() on exit? Default=FALSE
##' @param mapBubbles boolean (default=FALSE) add observation bubbles?
##' @param cutp optional vector of break points for the color scale on maps
##' @param map.pch pch for map points, default=16 (filled round).
##' @param ... Additional parameters for plot()
##' @return nothing
##' @export
##' @import maps mapdata tweedie
surveyIdxPlots<-function (x, dat, alt.idx = NULL, myids, cols = 1:length(x$pModels), 
    select = c("index", "map", "residuals", "fitVsRes"), par = list(mfrow = c(3, 
        3)), colors = rev(heat.colors(6)), map.cex = 1, plotByAge = TRUE, 
    legend = TRUE, predD = NULL, year = NULL, main=NULL, legend.signif=3,legend.pos="topright",restoreOldPar=FALSE,mapBubbles=FALSE,cutp = NULL,map.pch=16, ...) 
{
    if (!plotByAge & !is.null(par)){ 
        op<-par(par)
        if(restoreOldPar) on.exit(par(op))
    }
    mainwasnull <- is.null(main)
    for (a in cols) {
        if(mainwasnull) main <- paste("Age group", colnames(dat$Nage)[a])
        if (plotByAge & !is.null(par)){ 
            op<-par(par)
            if(restoreOldPar) on.exit(par(op))
        }
        if (any(select == "index")) {
            ys = range(as.numeric(levels(dat$Year)))
            ys = ys[1]:ys[2]
            yl = range(c(x$idx[, a],0,x$lo[,a],x$up[,a])/mean(x$idx[, a]),na.rm=TRUE) 
            if (!is.null(alt.idx) && a <= ncol(alt.idx)) {
                yl = range(c(alt.idx[, a]/mean(alt.idx[, a]), 
                  yl)) * 1.1
                plot(ys, alt.idx[, a]/mean(alt.idx[, a], na.rm = TRUE), 
                  ylim = yl, col = 2, ylab = "Index", xlab = "Year", 
                  main = main)
            }
            else {
                plot(ys, rep(NA, length(ys)), ylim = yl, col = 2, 
                  ylab = "Index", xlab = "Year", main = main )
            }
            idx = x$idx
            lo = x$lo
            up = x$up
            idx[x$idx <= 0] = NA
            lo[x$idx <= 0] = NA
            up[x$idx <= 0] = NA
            lines(ys, idx[, a]/mean(idx[, a], na.rm = TRUE), 
                lwd = 2)
            lines(ys, lo[, a]/mean(idx[, a], na.rm = TRUE), lwd = 2, 
                lty = 2)
            lines(ys, up[, a]/mean(idx[, a], na.rm = TRUE), lwd = 2, 
                lty = 2)
            if (legend && !is.null(alt.idx)) 
                legend(legend.pos, pch = c(1, NA), lty = c(NA, 
                  1), col = c(2, 1), legend = c("alt.idx", "GAM"))
        }
        if (any(select == "map")) {
            xlims = range(dat$lon, na.rm = TRUE)
            ylims = range(dat$lat, na.rm = TRUE)
            mapvals = NULL
            if (is.null(predD)) {
                tmp = subset(dat, haul.id %in% myids)
            }
            else {
                tmp = predD
            }
            if (is.null(year)) {
                concT = surveyIndex:::concTransform(log(x$gPreds[[a]]))
                mapvals = x$gPreds[[a]]
            }
            else {
                y = which(as.numeric(as.character(names(x$gPreds2[[a]]))) == 
                  year)
                if (length(y) == 0) 
                  stop(paste("Year", year, "age group", a, "not found."))
                concT = surveyIndex:::concTransform(log(x$gPreds2[[a]][[y]]))
                mapvals = x$gPreds2[[a]][[y]]
            }
            if (length(colors) > 1) 
                zFac = cut(concT, 0:length(colors)/length(colors))
            else zFac = 1
            if (length(map.cex) > 1) 
                sFac = cut(log(x$gPreds[[a]]), length(map.cex))
            else sFac = 1
            myCols = colors
            plot(tmp$lon, y = tmp$lat, col = 1, pch = 1, cex = map.cex[sFac], 
                xlim = xlims, ylim = ylims, xlab = "Longitude", 
                ylab = "Latitude", main = main, ...)
            points(tmp$lon, y = tmp$lat, col = myCols[zFac], 
                pch = 16, cex = map.cex[sFac])
            maps::map("worldHires", xlim = xlims, ylim = ylims, 
                fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
            if (legend){
                maxcuts = aggregate(mapvals ~ zFac, FUN=max)
                mincuts = aggregate(mapvals ~ zFac, FUN=min)
                mm = mean(mapvals)
                ml = signif(mincuts[,2]/mm,legend.signif)
                ml[1] = 0
                leg = paste0("[",ml,",",signif(maxcuts[,2]/mm,legend.signif),"]")
                legend(legend.pos, legend = leg, pch = 16, col = colors, bg = "white")
            }
        }
        if (any(select == "absolutemap") || any(select == "CVmap")) {
            if(is.null(year) || length(year)==0) stop("argument 'year' must be vector of length>=1 for type 'absolutemap'")
            if( !all(year %in% levels(dat$Year)) ) stop("invalid years selected")
            
            xlims = range(dat$lon, na.rm = TRUE)
            ylims = range(dat$lat, na.rm = TRUE)

            if(any(select == "absolutemap")) colsel = "gPreds2" else colsel = "gPreds2.CV" 
            goodyears = intersect(year,names(x[[colsel]][[a]])) 
            ## collect all years as data.frame
            ally = data.frame(val = numeric(0), year = character(0))
            for(y in goodyears){
                ally = rbind(ally, data.frame(val=x[[colsel]][[a]][[y]],
                                              year=y))
            }
            ally$conc = surveyIndex:::concTransform(log(ally$val))
            if(is.null(cutp)){
                ally$zFac=cut( ally$conc,0:length(colors)/length(colors))
            } else {
                if(length(cutp) != length(colors) + 1) stop("incompatible number of colors and length of cutp") 
                ally$zFac=cut( ally$val,cutp)
            }
            bubbleScale = 0.005*max(dat$Nage[,a])
            for(yy in year){
                sel = which(ally$year == yy)
                if (is.null(predD)) {
                    tmp = subset(dat, haul.id %in% myids)
                }
                else {
                    tmp = predD
                    if(is.list(tmp) && !class(tmp)%in%c("data.frame","DATRASraw")) tmp = predD[[as.character(yy)]]
                }
                
                plot(tmp$lon,y=tmp$lat,col=1,pch=1,cex=map.cex,xlab="Longitude",ylab="Latitude",axes=FALSE,type=ifelse(length(sel)>0,"p","n"))
                box()
                title(yy,line=1)
                if(length(sel)==0) next;
                points(tmp$lon,y=tmp$lat,col=colors[as.numeric(ally$zFac[sel])],pch=map.pch,cex=map.cex)
                maps::map('worldHires',xlim=xlims,ylim=ylims,fill=TRUE,plot=TRUE,add=TRUE,col=grey(0.5))
                if(mapBubbles){
                    dy = subset(dat,Year==yy)
                    points(dy$lon,dy$lat,cex=sqrt(dy$Nage[,a]/bubbleScale))
                }
                if (legend && yy==year[1]){
                    if(is.null(cutp)){
                        maxcuts = aggregate(val ~ zFac, data=ally, FUN=max)
                        mincuts = aggregate(val ~ zFac, data=ally, FUN=min)
                        mm = mean(ally$val)
                        ml = signif(mincuts[,2]/mm,legend.signif)
                        ml[1] = 0
                        leg = paste0("[",ml,",",signif(maxcuts[,2]/mm,legend.signif),"]")
                        legend(legend.pos, legend = leg, pch = 16, col = colors, bg = "white")
                    } else {
                        legend(legend.pos, legend = levels(ally$zFac), pch = 16, col = colors, bg = "white")          }
                }
            }##rof year
        }## absolutemap
        
        for (k in 1:length(select)) {
            ss = suppressWarnings(as.numeric(select[k]))
            if (!is.na(ss)) {
                plot.gam(x$pModels[[a]], select = ss, main = main, ...)
            }
        }
        if (any(select == "residuals") || any(select == "fitVsRes") || 
            any(select == "resVsYear") || any(select == "resVsShip") || 
            any(select == "spatialResiduals")) {
            
            resi <- x$residuals[[a]] ##residuals(x,a)
        }
        if (any(select == "residuals")) {
            hist(resi, nclass = 30, main = main, xlab = "Residuals",...)
        }
        if (any(select == "fitVsRes")) {
            plot(fitted(x$pModels[[a]]), residuals(x$pModels[[a]]), xlab = "Fitted", 
                ylab = "Residuals", main = main,...) ## TODO: use quantile residuals here too
        }
        if (any(select == "resVsYear")) {
            plot(dat$Year, resi, main = main, xlab = "Year", ylab = "Residuals", 
                ...)
        }
        if (any(select == "resVsShip")) {
            plot(dat$Ship, resi, main = main, xlab = "Year", ylab = "Residuals", 
                ...)
        }
        if (any(select == "spatialResiduals")) {
            scale <- 3 * map.cex
            if (is.null(year) || length(year)>1) 
                stop("a single year must be supplied")
            sel <- which(dat[[2]]$Year == as.character(year))
            plot(dat$lon, dat$lat, type = "n", 
                xlab = "Longitude", ylab = "Latitude", main = main,...)
            maps::map("worldHires", fill = TRUE, plot = TRUE, 
                add = TRUE, col = grey(0.5))
            positive = resi[sel] > 0
            points(dat$lon[sel][positive], dat$lat[sel][positive], 
                pch = 1, cex = scale * sqrt(resi[sel][positive]), 
                col = "blue")
            points(dat$lon[sel][!positive], dat$lat[sel][!positive], 
                pch = 1, cex = scale * sqrt(-resi[sel][!positive]), 
                col = "red")
        }
    }
}

##' Plot values for each color level in a map as created by 'surveyIdxPlots'.
##'
##' @title Plot values for each color level in a map as created by 'surveyIdxPlots'.
##' @param x object of class 'surveyIdx'
##' @param year vector of years (should be the same as used for 'surveyIdxPlots')
##' @param colors vector of colors to use
##' @param age age group (integer >=1)
##' @param legend.signif number of significant digits shown
##' @param cutp optional custom vector of cut points for color scales.
##' @param lwd line width for bars
##' @param las orientation of x-axis lables
##' @export
mapLegend<-function(x,year,colors = rev(heat.colors(6)),age = 1,legend.signif=3,cutp = NULL,lwd=10,las=2){
    colsel = "gPreds2"
    a = age
    goodyears = intersect(year,names(x[[colsel]][[a]])) 
    
    ally = data.frame(val = numeric(0), year = character(0))
    for(y in goodyears){
        ally = rbind(ally, data.frame(val=x[[colsel]][[a]][[y]],year=y))
    }
    ally$conc = surveyIndex:::concTransform(log(ally$val))
    if(is.null(cutp)){
        ally$zFac=cut( ally$conc,0:length(colors)/length(colors))
    } else {
        if(length(cutp) != length(colors) + 1) stop("incompatible number of colors and length of cutp") 
        ally$zFac=cut( ally$val,cutp)
    }
    mvals = aggregate(val ~ zFac, data=ally, FUN=mean)
    mm = mean(ally$val)
    if(is.null(cutp)){
        maxcuts = aggregate(val ~ zFac, data=ally, FUN=max)
        mincuts = aggregate(val ~ zFac, data=ally, FUN=min)
        ml = signif(mincuts[,2]/mm,legend.signif)
        ml[1] = 0
        leg = paste0("[",ml,",",signif(maxcuts[,2]/mm,legend.signif),"]")
        plot(mvals[,2],type="h",lwd=lwd+lwd*0.1,col=1,axes=FALSE)
        points(mvals[,2],type="h",lwd=lwd,col=colors)
        axis(1,at=1:length(colors),labels=leg,las=las)
    } else {
        leg = levels(ally$zFac)
        plot(mvals[,2],type="h",lwd=11,col=1,axes=FALSE)
        points(mvals[,2],type="h",lwd=10,col=colors)
        axis(1,at=1:length(colors),labels=leg,las=las)
    }
    
}


##' Plot estimates of factor levels (or random effects) from a GAM model.
##' 
##' @param x A model of class 'gam' or 'glm'
##' @param name name of the factor
##' @param invlink inverse link function
##' @param levs optional custom factor level names 
##' @param ... extra arguments to 'plot'
##' @export
factorplot<-function(x, name, invlink=exp, levs=NULL,ylim=NULL,... ){
    sel<- grep( name, names(coef(x))) 
    est <- coef(x)[ sel ]
    sds <- sqrt(diag( vcov(x) ))[ sel ]
    lo <- invlink(est - 2*sds)
    hi <- invlink(est + 2*sds)
    if(is.null(ylim)) ylim <- range(c(lo,hi))
    xs <- 1:length(est)

    nams <- names(est)
    if(!is.null(levs)) nams <- levs
    op <- par(las=2)
    on.exit(par(op))
    plot( xs, invlink(est), ylim=ylim,...,xaxt="n",xlab="", ylab="Estimate")
    axis(1,labels=nams,at=xs)
    arrows( xs,lo,y1=hi,angle=90,code=3,length=0.1)
}

##' Calculate and plot the distribution as a function of depth
##' 
##' @title Calculate and plot the distribution as a function of depth
##' @param m object of class 'surveyIdx'
##' @param grid grid object used when fitting 'm'
##' @param by depth bin size
##' @param type character vector indicating what to plot. 'dist' or 'mean' or both.
##' @param colfun color function to use for plotting.
##' @param lwd plotting line width
##' @param rangeFrac Fraction ]0;1[ indicating the percentile range in the mean plot. Default is 0.5, i.e. interquartile range.  
##' @param age which age group to plot
##' @param ... extra arguments to 'plot
##' @return a matrix with depth distribution in each year.
##' @export
depthDist<-function(m,grid,by=5,type=c("dist","mean"),colfun=colorRampPalette(c("red","lightgrey","blue")),lwd=1.5,rangeFrac=0.5,age=1,...){
    br = seq(min(grid$Depth)-1,max(grid$Depth)+1,by=by)
    gridd = cut(grid$Depth,breaks=br)
    ##dds = lapply(m$gPreds2,function(x) aggregate(x,by=list(gridd),FUN=mean))
    dds = lapply(m$gPreds2[[age]],function(x) xtabs(x ~ gridd)/sum(x))
    ddsrs = do.call(cbind,dds)
    
    if("dist" %in% type){
        cols = colfun(ncol(ddsrs))
        matplot(ddsrs,type="l",lty=1,col=cols,axes=FALSE,xlab="Depth",ylab="Proportion",lwd=lwd,...)
        axis(1,labels=levels(gridd),at=1:nrow(ddsrs))
        axis(2)
        box()
        sel = seq(1,nrow(m$idx),length.out=min(nrow(m$idx),5))
        legend("topright",lty=1,legend=rownames(m$idx)[sel],col=cols[sel],lwd=lwd)
    }
    if("mean" %in% type){
        md = br[-1]+0.5 ##aggregate(grid$Depth~gridd,FUN=mean)
        mdy = colSums(ddsrs*md)
        lo = md[apply( ddsrs,2,function(x) match(TRUE,cumsum(x)>rangeFrac/2) )]
        hi = md[apply( ddsrs,2,function(x) match(TRUE,cumsum(x)>(1-rangeFrac/2) ))]
        plot(rownames(m$idx),mdy,type="p",pch=16,cex=1.5,ylab="Mean depth",xlab="Year",ylim=range(c(lo,hi)),col=cols)
        arrows(as.numeric(rownames(m$idx)),y0=lo,y1=hi,code=3,angle=90,length=0.1,col=cols,lwd=lwd)
        abline(h=mean(mdy),lty=2)
    }
    ddsrs
}
