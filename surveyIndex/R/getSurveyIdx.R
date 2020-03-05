##' Calculate survey indices by age.
##'
##' This is based on the methods described in
##' Berg et al. (2014): "Evaluation of alternative age-based methods for estimating relative abundance from survey data in relation to assessment models",
##' Fisheries Research 151(2014) 91-99.
##' @title Calculate survey indices by age.
##' @param x DATRASraw object
##' @param ages vector of ages
##' @param myids haul.ids for grid
##' @param kvecP vector with spatial smoother max. basis dimension for each age group, strictly positive part of model 
##' @param kvecZ vector with spatial smoother max. basis dimension for each age group, presence/absence part of model (ignored for Tweedie models) 
##' @param gamma model degress of freedom inflation factor (see 'gamma' argument to gam() ) 
##' @param cutOff treat observations below this value as zero
##' @param fam distribution, either "Gamma","LogNormal", or "Tweedie".
##' @param useBIC use BIC for smoothness selection (overrides 'gamma' argument)
##' @param nBoot number of bootstrap samples used for calculating index confidence intervals
##' @param mc.cores number of cores for parallel processing
##' @param method smoothness selection method used by 'gam'
##' @param predD optional DATRASraw object, defaults to NULL. If not null this is used as grid.
##' @param modelZ vector of model formulae for presence/absence part, one pr. age group (ignored for Tweedie models)
##' @param modelP vector of model formulae for strictly positive repsonses, one pr. age group
##' @param knotsP optional list of knots to gam, strictly positive repsonses
##' @param knotsZ optional list of knots to gam, presence/absence
##' @param predfix optional named list of extra variables (besides Gear, HaulDur, Ship, and TimeShotHour),  that should be fixed during prediction step (standardized)
##' @param linkZ link function for the binomial part of the model, default: "logit" (not used for Tweedie models).
##' @param ... Optional extra arguments to "gam"
##' @return A survey index (list)
##' @examples
##' \dontrun{
##' library(surveyIndex)
##' ##downloadExchange("NS-IBTS",1994:2014)
##' dAll<-readExchangeDir(".",strict=FALSE)
##' mc.cores<-2; library(parallel)
##' d<-subset(dAll, Species=="Pollachius virens",Quarter==1,HaulVal=="V",StdSpecRecCode==1, Gear=="GOV")
##' dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
##' d<-addSpectrum(d,by=1)
##' ## get idea about number of age groups to include
##' agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
##' agetab.df<-as.data.frame(agetab)
##' ages<-1:8
##' ## require at least 1 aged individual in each year
##' for(a in ages){
##'     if(any(agetab.df$Freq[agetab.df$Age==a]<1))
##'         d<-fixAgeGroup(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
##' }
##' d<-subset(d,Age>=min(ages))
##' 
##' ###############################
##' ## Convert to numbers-at-age
##' ###############################
##' d.ysplit <- split(d, d$Year)
##' ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
##'                varCof=FALSE,maxK=50,mc.cores=mc.cores)
##' Nage<-mclapply(ALK,predict,mc.cores=mc.cores)
##' for(i in 1:length(ALK)) d.ysplit[[i]]$Nage=Nage[[i]];
##' dd <- do.call("c",d.ysplit)
##' 
##' ##############
##' ## Fit model
##' ##############
##' grid <- getGrid(dd, nLon=40)
##' ## set max basis dim for spatial smooths by age, P=positive and Z=zero/absence.
##' ## These are set relatively low here to speed up the example
##' kvP <- c(50,50,50,40,30,rep(10,length(ages)-5))
##' kvZ <- kvP / 2;
##' mP <- rep("Year+s(lon,lat,k=kvecP[a],bs='ts')+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(ages)  );
##' mZ <- rep("Year+s(lon,lat,k=kvecZ[a],bs='ts')+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(ages)  );
##' 
##' SIQ1 <- getSurveyIdx(dd,ages=ages,myids=grid[[3]],cutOff=0.1,kvecP=kvP,kvecZ=kvZ,
##'          modelZ=mZ,modelP=mP,mc.cores=mc.cores) ## if errors are encountered, debug with mc.cores=1 
##' 
##' strat.mean<-getSurveyIdxStratMean(dd,ages)
##' 
##' ## plot indices, distribution map, and estimated depth effects
##' surveyIdxPlots(SIQ1,dd,cols=ages,alt.idx=strat.mean,grid[[3]],par=list(mfrow=c(3,3)),legend=FALSE,
##'                select="index",plotByAge=FALSE)
##' 
##' surveyIdxPlots(SIQ1,dd,cols=ages,alt.idx=NULL,grid[[3]],par=list(mfrow=c(3,3)),legend=FALSE,
##'                 colors=rev(heat.colors(8)),select="map",plotByAge=FALSE)
##' 
##' surveyIdxPlots(SIQ1,dd,cols=ages,alt.idx=NULL,grid[[3]],par=list(mfrow=c(3,3)),
##'                 legend=FALSE,select="2",plotByAge=FALSE)
##' 
##' ## Calculate internal concistency and export to file
##' internalCons(SIQ1$idx)
##' exportSI(SIQ1$idx,ages=ages,years=levels(dd$Year),toy=mean(dd$timeOfYear),file="out.dat",
##'          nam="Survey index demo example")
##' }
##' @importFrom MASS mvrnorm
##' @export
getSurveyIdx <-
    function(x,ages,myids,kvecP=rep(12*12,length(ages)),kvecZ=rep(8*8,length(ages)),gamma=1.4,cutOff=1,fam="Gamma",useBIC=FALSE,nBoot=1000,mc.cores=1,method="ML",predD=NULL,
             modelZ=rep("Year+s(lon,lat,k=kvecZ[a],bs='ts')+s(Ship,bs='re',by=dum)+s(Depth,bs='ts')+s(TimeShotHour,bs='cc')",length(ages)  ),modelP=rep("Year+s(lon,lat,k=kvecP[a],bs='ts')+s(Ship,bs='re',by=dum)+s(Depth,bs='ts')+s(TimeShotHour,bs='cc')",length(ages)  ),knotsP=NULL,knotsZ=NULL,predfix=NULL,linkZ="logit", ...
             ){
        
        if(is.null(x$Nage)) stop("No age matrix 'Nage' found.");
        if(is.null(colnames(x$Nage))) stop("No colnames found on 'Nage' matrix.");
        if(length(modelP)<length(ages)) stop(" length(modelP) < length(ages)");
        if(length(kvecP)<length(ages)) stop(" length(kvecP) < length(ages)");
        if(fam[1]!="Tweedie"){ 
            if(length(modelZ)<length(ages)) stop(" length(modelZ) < length(ages)");
            if(length(kvecZ)<length(ages)) stop(" length(kvecZ) < length(ages)");
        }
        ## check for valid family names
        stopifnot(fam[1] %in% c("Gamma","LogNormal","Tweedie","negbin"))
        if(length(fam)<length(ages)) { famVec = rep(fam[1],length(ages)) } else famVec=fam;

        dataAges <- as.numeric(gsub("[+]","",colnames(x$Nage)))
        if(!all(ages%in%dataAges)) stop(paste0("age(s) ",setdiff(ages,dataAges)," not found in 'Nage' matrix"));
        x[[1]]$Year=as.factor(x[[1]]$Year);
        x[[2]]$Year=as.factor(x[[2]]$Year);
        pModels=list()
        zModels=list()
        gPreds=list() ##last data year's predictions
        gPreds2=list() ## all years predictions
        allobs=list() ## response vector (zeroes and positive)
        
        yearNum=as.numeric(as.character(x$Year));
        yearRange=min(yearNum):max(yearNum);

        ## Choose most frequent gear as reference gear
        gearNames=names(xtabs(~Gear,data=x[[2]]))
        myGear=names(xtabs(~Gear,data=x[[2]]))[which.max(xtabs(~Gear,data=x[[2]]))]
        
        resMat=matrix(NA,nrow=length(yearRange),ncol=length(ages));
        upMat=resMat;
        loMat=resMat;
        do.one.a<-function(a){
            age = which(dataAges==ages[a])
            ddd=x[[2]]; ddd$dum=1.0;
            ddd$A1=ddd$Nage[,age]
            gammaPos=gamma;
            gammaZ=gamma;
            if(useBIC){
                nZ=nrow(ddd);
                nPos=nrow(subset(ddd,A1>cutOff));
                gammaPos=log(nPos)/2;
                gammaZ=log(nZ)/2;
                cat("gammaPos: ",gammaPos," gammaZ: ",gammaZ,"\n");
            }
            pd = subset(ddd,A1>cutOff)
            if(famVec[a]=="LogNormal"){
                f.pos = as.formula( paste( "log(A1) ~",modelP[a]));
                f.0 = as.formula( paste( "A1>",cutOff," ~",modelZ[a]));
                
                print(system.time(m.pos<-tryCatch.W.E(gam(f.pos,data=subset(ddd,A1>cutOff),gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));

                if(class(m.pos)[2] == "error") {
                    print(m.pos)
                    stop("Error occured for age ", a, " in the positive part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
                }
                
                print(system.time(m0<-tryCatch.W.E(gam(f.0,gamma=gammaZ,data=ddd,family=binomial(link=linkZ),method=method,knots=knotsZ,na.action=na.fail,...))$value));

                if(class(m0)[2] == "error") {
                    print(m0)
                    stop("Error occured for age ", a, " in the binomial part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
                }
                
            } else if(famVec[a]=="Gamma"){
                f.pos = as.formula( paste( "A1 ~",modelP[a]));
                f.0 = as.formula( paste( "A1>",cutOff," ~",modelZ[a]));
                
                print(system.time(m.pos<-tryCatch.W.E(gam(f.pos,data=subset(ddd,A1>cutOff),family=Gamma(link="log"),gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));

                if(class(m.pos)[2] == "error") {
                    print(m.pos)
                    stop("Error occured for age ", a, " in the positive part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
                }
                
                print(system.time(m0<-tryCatch.W.E(gam(f.0,gamma=gammaZ,data=ddd,family=binomial(link=linkZ),method=method,knots=knotsZ,na.action=na.fail,...))$value));

                if(class(m0)[2] == "error") {
                    print(m0)
                    stop("Error occured for age ", a, " in the binomial part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
                }
            } else if(famVec[a]=="Tweedie"){
                ddd$A1[ ddd$A1<cutOff ] = 0
                pd = ddd
                f.pos = as.formula( paste( "A1 ~",modelP[a]));
                print(system.time(m.pos<-tryCatch.W.E(gam(f.pos,data=ddd,family=tw,gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));
                if(class(m.pos)[2] == "error") {
                    print(m.pos)
                    stop("Error occured for age ", a, ".\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
                }
                m0=NULL;
            } else if(famVec[a]=="negbin"){
                pd = ddd
                f.pos = as.formula( paste( "A1 ~",modelP[a]));
                print(system.time(m.pos<-tryCatch.W.E(gam(f.pos,data=ddd,family=nb,gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));
                if(class(m.pos)[2] == "error") {
                    print(m.pos)
                    stop("Error occured for age ", a, ".\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
                }
                m0=NULL;
            }
            ## Calculate total log-likelihood
            if(famVec[a]=="Tweedie" || famVec[a]=="negbin"){
                totll = logLik(m.pos)[1]
            } else {
                p0p =(1-predict(m0,type="response"))
                ppos=p0p[ddd$A1>cutOff] 
                p0m1=p0p[ddd$A1<=cutOff]
                if(famVec[a]=="Gamma")  totll=sum(log(p0m1))+sum(log(1-ppos))+logLik(m.pos)[1];
                ## if logNormal model, we must transform til log-likelihood to be able to use AIC
                ## L(y) = prod( dnorm( log y_i, mu_i, sigma^2) * ( 1 / y_i ) ) => logLik(y) = sum( log[dnorm(log y_i, mu_i, sigma^2)]  - log( y_i ) )
                if(famVec[a]=="LogNormal") totll=sum(log(p0m1))+ sum(log(1-ppos)) + logLik(m.pos)[1] - sum(m.pos$y);
            }
            
            if(is.null(predD)) predD=subset(ddd,haul.id %in% myids);
            res=numeric(length(yearRange));
            lores=res;
            upres=res;
            gp2=list()
            
            for(y in levels(ddd$Year)){ 
                ## take care of years with all zeroes
                if(!any(ddd$A1[ddd$Year==y]>cutOff)){
                    res[which(as.character(yearRange)==y)]=0;
                    upres[which(as.character(yearRange)==y)] = 0;
                    lores[which(as.character(yearRange)==y)] = 0;
                    next;
                }

                ## OBS: effects that should be removed should be included here
                predD$Year=y; predD$dum=0;
                predD$ctime=as.numeric(as.character(y));
                predD$TimeShotHour=mean(ddd$TimeShotHour)
                predD$Ship=names(which.max(summary(ddd$Ship)))
                predD$timeOfYear=mean(ddd$timeOfYear);
                predD$HaulDur=30.0
                predD$Gear=myGear;
                if(!is.null(predfix)){ ##optional extra variables for standardization
                    stopifnot(is.list(predfix))
                    for(n in names(predfix)){
                        predD[,n] = predfix[[n]]
                    }
                }

                p.1=try(predict(m.pos,newdata=predD,newdata.guaranteed=TRUE));
                if(!famVec[a] %in% c("Tweedie","negbin")) p.0=try(predict(m0,newdata=predD,type="response",newdata.guaranteed=TRUE));
                ## take care of failing predictions
                if(!is.numeric(p.1) | (!famVec[a] %in% c("Tweedie","negbin") && !is.numeric(p.0))) {
                    res[which(as.character(yearRange)==y)]=0;
                    upres[which(as.character(yearRange)==y)] = 0;
                    lores[which(as.character(yearRange)==y)] = 0;
                    next;
                }
                sig2=m.pos$sig2;
                
                if(famVec[a]=="Gamma") { res[which(as.character(yearRange)==y)] = sum(p.0*exp(p.1)); gPred=p.0*exp(p.1) }
                if(famVec[a]=="LogNormal")  { res[which(as.character(yearRange)==y)] = sum(p.0*exp(p.1+sig2/2)); gPred=p.0*exp(p.1+sig2/2) }
                if(famVec[a] %in% c("Tweedie","negbin"))  { res[which(as.character(yearRange)==y)] = sum(exp(p.1)); gPred=exp(p.1) }
                gp2[[y]]=gPred;
                if(nBoot>10){
                    Xp.1=predict(m.pos,newdata=predD,type="lpmatrix");
                    brp.1=mvrnorm(n=nBoot,coef(m.pos),m.pos$Vp);
                    if(!famVec[a] %in% c("Tweedie","negbin")){
                        Xp.0=predict(m0,newdata=predD,type="lpmatrix");
                        brp.0=mvrnorm(n=nBoot,coef(m0),m0$Vp);
                        OS0 = matrix(0,nrow(predD),nBoot);
                        terms.0=terms(m0)
                        if(!is.null(m0$offset)){
                            off.num.0 <- attr(terms.0, "offset")
                            
                            for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0, 
                                                                        "variables")[[i + 1]], predD)
                        }
                        rep0=m0$family$linkinv(Xp.0%*%t(brp.0)+OS0);
                    }
                    OS.pos = matrix(0,nrow(predD),nBoot);
                    terms.pos=terms(m.pos)
                    if(!is.null(m.pos$offset)){
                        off.num.pos <- attr(terms.pos, "offset")
                        
                        for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos, 
                                                                            "variables")[[i + 1]], predD)
                    }
                    
                    
                    if(famVec[a]=="LogNormal"){
                        rep1=exp(Xp.1%*%t(brp.1)+sig2/2+OS.pos);
                    } else {
                        rep1=exp(Xp.1%*%t(brp.1)+OS.pos);
                    } 

                    if(!famVec[a] %in% c("Tweedie","negbin")){
                        idxSamp = colSums(rep0*rep1);
                    } else {
                        idxSamp = colSums(rep1);
                    }
                    upres[which(as.character(yearRange)==y)] = quantile(idxSamp,0.975);
                    lores[which(as.character(yearRange)==y)] = quantile(idxSamp,0.025);
                }
            } ## rof years
            list(res=res,m.pos=m.pos,m0=m0,lo=lores,up=upres,gp=gPred,ll=totll,pd=pd,gp2=gp2);
        }## end do.one
        noAges=length(ages);
        rr=parallel::mclapply(1:noAges,do.one.a,mc.cores=mc.cores);
        logl=0;
        for(a in 1:noAges){
            resMat[,a]=rr[[a]]$res;
            zModels[[a]]=rr[[a]]$m0;
            pModels[[a]]=rr[[a]]$m.pos;
            loMat[,a]=rr[[a]]$lo;
            upMat[,a]=rr[[a]]$up;
            gPreds[[a]]=rr[[a]]$gp;
            logl=logl+rr[[a]]$ll
            gPreds2[[a]]=rr[[a]]$gp2
            allobs[[a]]=x[[2]]$Nage[,a]
        }
        getEdf<-function(m) sum(m$edf)
        totEdf=sum( unlist( lapply(zModels,getEdf))) + sum( unlist( lapply(pModels,getEdf)));
        rownames(resMat)<-yearRange
        colnames(resMat)<-ages
        out <- list(idx=resMat,zModels=zModels,pModels=pModels,lo=loMat,up=upMat,gPreds=gPreds,logLik=logl,edfs=totEdf,gPreds2=gPreds2,family=famVec, cutOff=cutOff, dataAges=dataAges, yearNum=yearNum, refGear=myGear, predfix = predfix, knotsP=knotsP, knotsZ=knotsZ, allobs=allobs);
        class(out) <- "surveyIdx"
        out
    }


##' Re-compute standardized survey indices for an alternative grid from a previous fitted "surveyIdx" model.
##' 
##' @title Re-compute standardized survey indices for an alternative grid from a previous fitted "surveyIdx" model.
##' @param x DATRASraw dataset
##' @param model object of class "surveyIdx" as created by "getSurveyIdx"
##' @param predD optional DATRASraw object, defaults to NULL. If not null this is used as grid.
##' @param myids haul.ids for grid
##' @param nBoot number of bootstrap samples used for calculating index confidence intervals
##' @param predfix optional named list of extra variables (besides Gear, HaulDur, Ship, and TimeShotHour),  that should be fixed during prediction step (standardized)
##' @param mc.cores mc.cores number of cores for parallel processing
##' @return An object of class "surveyIdx"
##' @importFrom MASS mvrnorm
##' @export
redoSurveyIndex<-function(x,model,predD=NULL,myids,nBoot=1000,predfix,mc.cores=1){        
    ages = as.numeric(colnames(model$idx))
    dataAges <- model$dataAges
    famVec = model$family
    cutOff = model$cutOff
    
    yearNum=model$yearNum
    yearRange=min(yearNum):max(yearNum);

    gPreds=list() ##last data year's predictions
    gPreds2=list() ## all years predictions
    
    myGear=model$refGear 
    
    resMat=matrix(NA,nrow=length(yearRange),ncol=length(ages));
    upMat=resMat;
    loMat=resMat;
    do.one.a<-function(a){
        age = which(dataAges==ages[a])
        ddd=x[[2]]; ddd$dum=1.0;
        ddd$A1=ddd$Nage[,age]
        m.pos = model$pModels[[a]]
        m0 = NULL
        if(!famVec[a] %in% c("Tweedie","negbin")) m0 = model$zModels[[a]]
        if(is.null(predD)) predD=subset(ddd,haul.id %in% myids);
        res=numeric(length(yearRange));
        lores=res;
        upres=res;
        gp2=list()
        do.one.y<-function(y){
            
            cat("Doing year ",y,"\n")
            ## take care of years with all zeroes
            if(!any(ddd$A1[ddd$Year==y]>cutOff)){
                return(list(res=0,upres=0,lores=0,gp2=NULL))
            }

            ## OBS: effects that should be removed should be included here
            predD$Year=y; predD$dum=0;
            predD$ctime=as.numeric(as.character(y));
            predD$TimeShotHour=mean(ddd$TimeShotHour)
            predD$Ship=names(which.max(summary(ddd$Ship)))
            predD$timeOfYear=mean(ddd$timeOfYear);
            predD$HaulDur=30.0
            predD$Gear=myGear;
            if(!is.null(predfix)){ ##optional extra variables for standardization
                stopifnot(is.list(predfix))
                for(n in names(predfix)){
                    predD[,n] = predfix[[n]]
                }
            }

            p.1=try(predict(m.pos,newdata=predD,newdata.guaranteed=TRUE));
            if(!famVec[a] %in% c("Tweedie","negbin")) p.0=try(predict(m0,newdata=predD,type="response",newdata.guaranteed=TRUE));
            ## take care of failing predictions
            if(!is.numeric(p.1) | (!famVec[a] %in% c("Tweedie","negbin") && !is.numeric(p.0))) {
                return(list(res=0,upres=0,lores=0,gp2=NULL))
            }
            sig2=m.pos$sig2;
            idx = NA
            if(famVec[a]=="Gamma") { idx <- sum(p.0*exp(p.1)); gPred=p.0*exp(p.1) }
            if(famVec[a]=="LogNormal")  { idx <- sum(p.0*exp(p.1+sig2/2)); gPred=p.0*exp(p.1+sig2/2) }
            if(famVec[a] %in% c("Tweedie","negbin"))  { idx <- sum(exp(p.1)); gPred=exp(p.1) }
            
            if(nBoot>10){
                Xp.1=predict(m.pos,newdata=predD,type="lpmatrix");
                brp.1=mvrnorm(n=nBoot,coef(m.pos),m.pos$Vp);
                if(!famVec[a] %in% c("Tweedie","negbin")){
                    Xp.0=predict(m0,newdata=predD,type="lpmatrix");
                    brp.0=mvrnorm(n=nBoot,coef(m0),m0$Vp);
                    OS0 = matrix(0,nrow(predD),nBoot);
                    terms.0=terms(m0)
                    if(!is.null(m0$offset)){
                        off.num.0 <- attr(terms.0, "offset")
                        
                        for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0, 
                                                                    "variables")[[i + 1]], predD)
                    }
                    rep0=m0$family$linkinv(Xp.0%*%t(brp.0)+OS0);
                }
                OS.pos = matrix(0,nrow(predD),nBoot);
                terms.pos=terms(m.pos)
                if(!is.null(m.pos$offset)){
                    off.num.pos <- attr(terms.pos, "offset")
                    
                    for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos, 
                                                                        "variables")[[i + 1]], predD)
                }
                
                
                if(famVec[a]=="LogNormal"){
                    rep1=exp(Xp.1%*%t(brp.1)+sig2/2+OS.pos);
                } else {
                    rep1=exp(Xp.1%*%t(brp.1)+OS.pos);
                } 

                if(!famVec[a] %in% c("Tweedie","negbin")){
                    idxSamp = colSums(rep0*rep1);
                } else {
                    idxSamp = colSums(rep1);
                }
                
                return(list(res=idx,upres=quantile(idxSamp,0.975),lores=quantile(idxSamp,0.025),gp2=gPred))
            }
        } ## rof years
        yres = parallel::mclapply(levels(ddd$Year),do.one.y,mc.cores=mc.cores)
        for(y in levels(ddd$Year)) {
            ii = which(as.character(yearRange)==y) 
            res[ii] = yres[[ii]]$res
            upres[ii] = yres[[ii]]$upres
            lores[ii] = yres[[ii]]$lores
            gp2[[y]] = yres[[ii]]$gp2            
        }
        list(res=res,m.pos=m.pos,m0=m0,lo=lores,up=upres,gp=tail(gp2,1),gp2=gp2);
    }## end do.one
    noAges=length(ages);
    rr=lapply(1:noAges,do.one.a);
    logl=0;
    for(a in 1:noAges){
        resMat[,a]=rr[[a]]$res;
        loMat[,a]=rr[[a]]$lo;
        upMat[,a]=rr[[a]]$up;
        gPreds[[a]]=rr[[a]]$gp;
        gPreds2[[a]]=rr[[a]]$gp2
    }
    rownames(resMat)<-yearRange
    colnames(resMat)<-ages
    out <- list(idx=resMat,zModels=model$zModels,pModels=model$pModels,lo=loMat,up=upMat,gPreds=gPreds,logLik=model$logLik,edfs=model$edfs,pData=model$pData,gPreds2=gPreds2,
                family=famVec, cutOff=cutOff, dataAges=dataAges, yearNum=yearNum, refGear=myGear, predfix = predfix, knotsP=model$knotsP, knotsZ=model$knotsZ);
    class(out) <- "surveyIdx"
    out
}

