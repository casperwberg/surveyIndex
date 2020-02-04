##' Visualize results from a survey index model fitted with getSurveyIdx().
##'
##' @title Visualize results from a survey index model fitted with getSurveyIdx().
##' @param x Survey index as produced by getSurveyIndex()
##' @param dat DATRASraw object
##' @param alt.idx optional matrix with alternative index
##' @param myids vector of haul ids that constitute the grid
##' @param cols which age columns to consider?
##' @param select character vector of chosen plots. Either one of "index","map","absolutemap","residuals","fitVsRes",""resVsYear","resVsShip","spatialResiduals", or a number. Numbers refer to smooths in the order they appear in the formula.
##' @param par 'par' settings for plotting (a named list).
##' @param colors colors for spatial effect.
##' @param map.cex size of grid points on maps
##' @param plotByAge boolean (default=TRUE). If true, par(par) is called for each age group.
##' @param legend boolean (default=TRUE). add legends to plot?
##' @param predD DATRASraw object with grid (optional). Overrides 'myids' if supplied.
##' @param year numeric scalar or vector (default=NULL). If 'select' equals 'map' a specific year can be chosen (only meaningful for time-varying spatial effects). If select equals 'absolutemap' year must be a vector. 
##' @param main optional main title (overrides default title)
##' @param legend.signif Number of significant digits in map legends
##' @param legend.pos Position of legend (e.g. "bottomleft") see ?legend
##' @param ... Additional parameters for plot()
##' @return nothing
##' @export
##' @import maps mapdata tweedie
surveyIdxPlots<-function (x, dat, alt.idx = NULL, myids, cols = 1:length(x$pModels), 
    select = c("index", "map", "residuals", "fitVsRes"), par = list(mfrow = c(3, 
        3)), colors = rev(heat.colors(6)), map.cex = 1, plotByAge = TRUE, 
    legend = TRUE, predD = NULL, year = NULL, main=NULL, legend.signif=3,legend.pos="topright",...) 
{
    if (!plotByAge & !is.null(par)) 
        par(par)
    mainwasnull <- is.null(main)
    for (a in cols) {
        if(mainwasnull) main <- paste("Age group", colnames(dat$Nage)[a])
        if (plotByAge & !is.null(par)) 
            par(par)
        if (any(select == "index")) {
            ys = range(as.numeric(levels(dat$Year)))
            ys = ys[1]:ys[2]
            yl = range(c(x$idx[, a]/mean(x$idx[, a]),0)) * 1.1
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
            if (legend) 
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
        if (any(select == "absolutemap")) {
            if(is.null(year) || length(year)<1) stop("argument 'year' must be vector of length>=1 for type 'absolutemap'")
            if( !all(year %in% levels(dat$Year)) ) stop("invalid years selected")
            xlims = range(dat$lon, na.rm = TRUE)
            ylims = range(dat$lat, na.rm = TRUE)
            if (is.null(predD)) {
                tmp = subset(dat, haul.id %in% myids)
            }
            else {
                tmp = predD
            }
            ## collect all years as data.frame
            ally = data.frame(val=x$gPreds2[[a]][[1]],year=as.character(levels(dat$Year)[1]))
            cc=0
            for(y in levels(dat$Year)){
                cc=cc+1
                ally = rbind(ally, data.frame(val=x$gPreds2[[a]][[cc]],
                                              year=as.character(levels(dat$Year)[cc])))
            }
            ally$conc = surveyIndex:::concTransform(log(ally$val))
            ally$zFac=cut( ally$conc,0:length(colors)/length(colors))
            for(yy in year){
                plot(tmp$lon,y=tmp$lat,col=1,pch=1,cex=map.cex,xlab="Longitude",ylab="Latitude",axes=FALSE)
                box()
                title(yy,line=1)
                sel = which(ally$year==yy)
                points(tmp$lon,y=tmp$lat,col=colors[as.numeric(ally$zFac[sel])],pch=16,cex=map.cex)
                maps::map('worldHires',xlim=xlims,ylim=ylims,fill=TRUE,plot=TRUE,add=TRUE,col=grey(0.5))
                if (legend && yy==year[1]){
                    maxcuts = aggregate(val ~ zFac, data=ally, FUN=max)
                    mincuts = aggregate(val ~ zFac, data=ally, FUN=min)
                    mm = mean(ally$val)
                    ml = signif(mincuts[,2]/mm,legend.signif)
                    ml[1] = 0
                    leg = paste0("[",ml,",",signif(maxcuts[,2]/mm,legend.signif),"]")
                    legend(legend.pos, legend = leg, pch = 16, col = colors, bg = "white")
                }
            }
        }
        for (k in 1:length(select)) {
            ss = suppressWarnings(as.numeric(select[k]))
            if (!is.na(ss)) {
                plot.gam(x$pModels[[a]], select = ss, main = main, ...)
            }
        }
        if (any(select == "residuals") || any(select == "fitVsRes") || 
            any(select == "resVsYear") || any(select == "resVsShip") || 
            any(select == "spatialResiduals")) {
            ## if (pmatch("Tweedie", x$pModels[[a]]$family$family, 
            ##     nomatch = -1) == 1) {
            ##     resi <- qres.tweedie(x$pModels[[a]])
            ## }
            ## else {
            ##     resi <- residuals(x$pModels[[a]])
            ## }
            resi <- residuals(x,a)
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
            scale <- 3
            if (is.null(year)) 
                stop("a year must be supplied")
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

