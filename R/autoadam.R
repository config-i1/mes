#' @param parallel If TRUE, the estimation of ADAM models is done in parallel (used in \code{auto.adam} only).
#' If the number is provided (e.g. \code{parallel=41}), then the specified number of cores is set up.
#' WARNING! Packages \code{foreach} and either \code{doMC} (Linux and Mac only)
#' or \code{doParallel} are needed in order to run the function in parallel.
#' @param fast If \code{TRUE}, then some of the orders of ARIMA are
#' skipped in the order selection. This is not advised for models with \code{lags} greater than 12.#'
#' @examples
#' ourModel <- auto.adam(rnorm(100,100,10), model="ZZN", lags=c(1,4), orders=list(ar=c(2,2),ma=c(2,2),select=TRUE))
#'
#' @rdname adam
#' @export
auto.adam <- function(y, model="ZXZ", lags=c(frequency(y)), orders=list(ar=c(0),i=c(0),ma=c(0),select=FALSE),
                      distribution=c("dnorm","dlogis","dlaplace","ds",
                                     "dlnorm","dllaplace","dls","dinvgauss"),
                      h=0, holdout=FALSE,
                      persistence=NULL, phi=NULL, initial=c("optimal","backcasting"), arma=NULL,
                      occurrence=c("none","auto","fixed","general","odds-ratio","inverse-odds-ratio","direct"),
                      ic=c("AICc","AIC","BIC","BICc"), bounds=c("usual","admissible","none"),
                      xreg=NULL, xregDo=c("use","select","adapt"),
                      silent=TRUE, parallel=FALSE, fast=TRUE, ...){
    # Copyright (C) 2020 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    #### modelDo, ic ####
    if(any(unlist(strsplit(model,""))=="C")){
        modelDo <- "combine";
    }
    else{
        modelDo <- "select";
    }
    nModels <- length(distribution);

    ic <- match.arg(ic,c("AICc","AIC","BIC","BICc"));
    ICFunction <- switch(ic,
                         "AIC"=AIC,
                         "AICc"=AICc,
                         "BIC"=BIC,
                         "BICc"=BICc);

    # The function checks the provided parameters of adam and/or oes
    ##### data #####
    # If this is simulated, extract the actuals
    if(is.adam.sim(y) || is.smooth.sim(y)){
        y <- y$data;
    }
    # If this is Mdata, use all the available stuff
    else if(inherits(y,"Mdata")){
        h <- y$h;
        holdout <- TRUE;
        lags <- frequency(y$x);
        y <- ts(c(y$x,y$xx),start=start(y$x),frequency=frequency(y$x));
    }

    # Extract index from the object in order to use it later
    ### tsibble has its own index function, so shit happens becaus of it...
    if(inherits(y,"tbl_ts")){
        yIndex <- y[[1]];
        if(any(duplicated(yIndex))){
            warning(paste0("You have duplicated time stamps in the variable ",responseName,
                           ". We will refactor this."),call.=FALSE);
            yIndex <- yIndex[1] + c(1:length(y[[1]])) * diff(tail(yIndex,2));
        }
    }
    else{
        yIndex <- try(time(y),silent=TRUE);
        # If we cannot extract time, do something
        if(inherits(yIndex,"try-error")){
            if(!is.null(dim(y))){
                yIndex <- as.POSIXct(rownames(y));
            }
            else{
                yIndex <- c(1:length(y));
            }
        }
    }
    yClasses <- class(y);

    # If this is something like a matrix
    if(!is.null(ncol(y))){
        # If we deal with data.table / tibble / data.frame, the syntax is different.
        # We don't want to import specific classes, so just use inherits()
        if(inherits(y,"tbl_ts")){
            # With tsibble we cannot extract explanatory variables easily...
            y <- y$value;
        }
        else if(inherits(y,"data.table") || inherits(y,"tbl") || inherits(y,"data.frame")){
            if(ncol(y)>1){
                xreg <- y[,-1];
            }
            y <- y[[1]];
        }
        else if(inherits(y,"zoo")){
            if(ncol(y)>1){
                xreg <- as.data.frame(y[,-1]);
            }
            y <- zoo(y[,1],order.by=time(y));
        }
        else{
            if(ncol(y)>1){
                xreg <- y[,-1];
            }
            y <- y[,1];
        }
    }

    #### Create logical, determining, what we are dealing with
    # ETS
    etsModel <- all(model!="NNN");
    # ARIMA + ARIMA select
    if(is.list(orders)){
        arimaModel <- any(c(orders$ar,orders$i,orders$ma)>0);
    }
    else{
        arimaModel <- any(orders>0);
    }
    #### Checks of provided parameters for ARIMA selection ####
    if(arimaModel && is.list(orders)){
        arimaModelSelect <- orders$select;
        arMax <- orders$ar;
        iMax <- orders$i;
        maMax <- orders$ma;

        if(is.null(arimaModelSelect)){
            arimaModelSelect <- FALSE;
        }

        if(any(c(arMax,iMax,maMax)<0)){
            stop("Funny guy! How am I gonna construct a model with negative order?",call.=FALSE);
        }

        # If there are no lags for the basic components, correct this.
        if(sum(lags==1)==0){
            lags <- c(1,lags);
        }

        # If there are zero lags, drop them
        if(any(lags==0)){
            arMax <- arMax[lags!=0];
            iMax <- iMax[lags!=0];
            maMax <- maMax[lags!=0];
            lags <- lags[lags!=0];
        }

        # Define maxorder and make all the values look similar (for the polynomials)
        maxorder <- max(length(arMax),length(iMax),length(maMax));
        if(length(arMax)!=maxorder){
            arMax <- c(arMax,rep(0,maxorder-length(arMax)));
        }
        if(length(iMax)!=maxorder){
            iMax <- c(iMax,rep(0,maxorder-length(iMax)));
        }
        if(length(maMax)!=maxorder){
            maMax <- c(maMax,rep(0,maxorder-length(maMax)));
        }

        # If zeroes are defined as orders for some lags, drop them.
        if(any((arMax + iMax + maMax)==0)){
            orders2leave <- (arMax + iMax + maMax)!=0;
            if(all(!orders2leave)){
                orders2leave <- lags==min(lags);
            }
            arMax <- arMax[orders2leave];
            iMax <- iMax[orders2leave];
            maMax <- maMax[orders2leave];
            lags <- lags[orders2leave];
        }

        # Get rid of duplicates in lags
        if(length(unique(lags))!=length(lags)){
            lagsNew <- unique(lags);
            arMaxNew <- iMaxNew <- maMaxNew <- lagsNew;
            for(i in 1:length(lagsNew)){
                arMaxNew[i] <- max(arMax[which(lags==lagsNew[i])],na.rm=TRUE);
                iMaxNew[i] <- max(iMax[which(lags==lagsNew[i])],na.rm=TRUE);
                maMaxNew[i] <- max(maMax[which(lags==lagsNew[i])],na.rm=TRUE);
            }
            arMax <- arMaxNew;
            iMax <- iMaxNew;
            maMax <- maMaxNew;
            lags <- lagsNew;
        }

        # Order things, so we would deal with the lowest level of seasonality first
        arMax <- arMax[order(lags,decreasing=FALSE)];
        iMax <- iMax[order(lags,decreasing=FALSE)];
        maMax <- maMax[order(lags,decreasing=FALSE)];
        lags <- sort(lags,decreasing=FALSE);
    }
    else{
        arMax <- iMax <- maMax <- NULL;
        arimaModelSelect <- FALSE;
    }
    # xreg - either as a separate variable or as a matrix for y
    xregModel <- !is.null(xreg) || (!is.null(dim(y)) && ncol(y>0));

    # Check the parallel parameter and set the number of cores
    if(is.numeric(parallel)){
        nCores <- parallel;
        parallel <- TRUE
    }
    else{
        nCores <- min(parallel::detectCores() - 1, nModels);
    }

    # If this is parallel, then load the required packages
    if(parallel){
        if(!requireNamespace("foreach", quietly = TRUE)){
            stop("In order to run the function in parallel, 'foreach' package must be installed.", call. = FALSE);
        }
        if(!requireNamespace("parallel", quietly = TRUE)){
            stop("In order to run the function in parallel, 'parallel' package must be installed.", call. = FALSE);
        }

        # Check the system and choose the package to use
        if(Sys.info()['sysname']=="Windows"){
            if(requireNamespace("doParallel", quietly = TRUE)){
                cat(paste0("Setting up ", nCores, " clusters using 'doParallel'..."));
                cat("\n");
                cluster <- parallel::makeCluster(nCores);
                doParallel::registerDoParallel(cluster);
            }
            else{
                stop("Sorry, but in order to run the function in parallel, you need 'doParallel' package.",
                     call. = FALSE);
            }
        }
        else{
            if(requireNamespace("doMC", quietly = TRUE)){
                doMC::registerDoMC(nCores);
                cluster <- NULL;
            }
            else if(requireNamespace("doParallel", quietly = TRUE)){
                cat(paste0("Setting up ", nCores, " clusters using 'doParallel'..."));
                cat("\n");
                cluster <- parallel::makeCluster(nCores);
                doParallel::registerDoParallel(cluster);
            }
            else{
                stop(paste0("Sorry, but in order to run the function in parallel, you need either ",
                            "'doMC' (prefered) or 'doParallel' package."),
                     call. = FALSE);
            }
        }
    }

    if(!silent){
        if(!parallel){
            cat("Evaluating models with different distributions... ");
        }
        else{
            cat(paste0("Working..."));
        }
    }

    #### The function that does the loop and returns a list of ETS(X) models ####
    adamReturner <- function(y, model, lags, orders,
                             distribution, h, holdout,
                             persistence, phi, initial, arma,
                             occurrence, ic, bounds,
                             xreg, xregDo, parallel,
                             arimaModelSelect, arMax, iMax, maMax, ...){
        # If we select ARIMA, don't do it in the first step
        if(arimaModelSelect){
            ordersToUse <- c(0,0,0);
        }
        else{
            ordersToUse <- orders;
        }

        if(!parallel){
            # Prepare the list of models
            selectedModels <- vector("list",length(distribution));
            for(i in 1:length(distribution)){
                if(!silent){
                    cat(paste0(distribution[i],", "));
                }
                selectedModels[[i]] <- adam(y=y, model=model, lags=lags, orders=ordersToUse,
                                            distribution=distribution[i],
                                            h=h, holdout=holdout,
                                            persistence=persistence, phi=phi, initial=initial, arma=arma,
                                            occurrence=occurrence, ic=ic, bounds=bounds,
                                            xreg=xreg, xregDo=xregDo, silent=TRUE, ...);

                if(arimaModelSelect){
                    selectedModels[[i]] <- arimaSelector(y=y, model=model,
                                                         lags=lags, arMax=arMax, iMax=iMax, maMax=maMax,
                                                         distribution=selectedModels[[i]]$distribution, h=h, holdout=holdout,
                                                         persistence=persistence, phi=phi, initial=initial,
                                                         occurrence=occurrence, ic=ic, bounds=bounds, fast=fast,
                                                         silent=silent, xreg=xreg, xregDo=xregDo,
                                                         testModelETS=selectedModels[[i]], ...)
                }
            }
        }
        else{
            selectedModels <- foreach::`%dopar%`(foreach::foreach(i=1:length(distribution)),{
                testModel <- adam(y=y, model=model, lags=lags, orders=ordersToUse,
                                  distribution=distribution[i],
                                  h=h, holdout=holdout,
                                  persistence=persistence, phi=phi, initial=initial, arma=arma,
                                  occurrence=occurrence, ic=ic, bounds=bounds,
                                  xreg=xreg, xregDo=xregDo, silent=TRUE, ...)

                if(arimaModelSelect){
                    testModel <- arimaSelector(y=y, model=model,
                                               lags=lags, arMax=arMax, iMax=iMax, maMax=maMax,
                                               distribution=testModel$distribution, h=h, holdout=holdout,
                                               persistence=persistence, phi=phi, initial=initial,
                                               occurrence=occurrence, ic=ic,
                                               bounds=bounds, fast=fast,
                                               silent=TRUE, xreg=xreg, xregDo=xregDo,
                                               testModelETS=testModel, ...)
                }
                return(testModel);
            })
        }
        return(selectedModels);
    }

    #### ARIMA selection script ####
    if(arimaModelSelect){
        if(!is.null(arma)){
            warning("ARIMA order selection cannot be done with the provided arma parameters. Dropping them.",
                    call.=FALSE);
            arma <- NULL;
        }

        #### Function corrects IC taking number of parameters on previous step ####
        icCorrector <- function(icValue, ic, nParam, obsNonzero, nParamNew){
            if(ic=="AIC"){
                llikelihood <- (2*nParam - icValue)/2;
                correction <- 2*nParamNew - 2*llikelihood;
            }
            else if(ic=="AICc"){
                llikelihood <- (2*nParam*obsNonzero/(obsNonzero-nParam-1) - icValue)/2;
                correction <- 2*nParamNew*obsNonzero/(obsNonzero-nParamNew-1) - 2*llikelihood;
            }
            else if(ic=="BIC"){
                llikelihood <- (nParam*log(obsNonzero) - icValue)/2;
                correction <- nParamNew*log(obsNonzero) - 2*llikelihood;
            }
            else if(ic=="BICc"){
                llikelihood <- ((nParam*log(obsNonzero)*obsNonzero)/(obsNonzero-nParam-1) - icValue)/2;
                correction <- (nParamNew*log(obsNonzero)*obsNonzero)/(obsNonzero-nParamNew-1) - 2*llikelihood;
            }

            return(correction);
        }

        #### The function that selects ARIMA orders for the provided data ####
        arimaSelector <- function(y, model, lags, arMax, iMax, maMax,
                                  distribution, h, holdout,
                                  persistence, phi, initial,
                                  occurrence, ic, bounds, fast,
                                  silent, xreg, xregDo, testModelETS, ...){
            silentDebug <- TRUE;

            # Save the original values
            modelOriginal <- model;
            xregOriginal <- xreg;
            occurrenceOriginal <- occurrence;
            persistenceOriginal <- persistence;
            phiOriginal <- phi;

            # If the ETS model was done before this, then extract residuals
            if(is.adam(testModelETS)){
                dataAR <- dataI <- dataMA <- yInSample <- residuals(testModelETS);
                model <- "NNN";
                xreg <- NULL;
                occurrence <- "none"
                persistence <- NULL;
                phi <- NULL;

                # Don't count the scale term
                nParamOriginal <- nparam(testModelETS)-1;
            }
            else{
                # Fit Naive and get the parameters
                testModelETS <- adam(y,model="ANN",lags=1,distribution=distribution,
                                     h=h,holdout=holdout,persistence=0,initial=mean(y),
                                     occurrence=occurrence,bounds="none",silent=TRUE);
                dataAR <- dataI <- dataMA <- yInSample <- actuals(testModelETS);

                # Originally, we only have a constant
                nParamOriginal <- 1;
            }
            testModel <- testModelETS;
            bestIC <- bestICI <- ICFunction(testModel);
            obsNonzero <- nobs(testModelETS,all=FALSE);

            if(silentDebug){
                cat("Best IC: "); cat(bestIC); cat("\n");
            }
            if(!silent){
                cat(" Selecting ARIMA orders...    ");
            }

            # 1 stands for constant/no constant, another one stands for ARIMA(0,0,0)
            if(all(maMax==0)){
                nModelsARIMA <- prod(iMax + 1) * (1 + sum(arMax));
            }
            else{
                nModelsARIMA <- prod(iMax + 1) * (1 + sum(maMax*(1 + sum(arMax))));
            }
            ICValue <- 1E+100;
            m <- 0;

            lagsTest <- maTest <- arTest <- rep(0,length(lags));
            arBest <- maBest <- iBest <- rep(0,length(lags));
            arBestLocal <- maBestLocal <- arBest;

            iOrders <- matrix(0,prod(iMax+1),ncol=length(iMax));

            ##### Loop for differences #####
            # Prepare table with differences
            if(any(iMax!=0)){
                iOrders[,1] <- rep(c(0:iMax[1]),times=prod(iMax[-1]+1));
                if(length(iMax)>1){
                    for(seasLag in 2:length(iMax)){
                        iOrders[,seasLag] <- rep(c(0:iMax[seasLag]),each=prod(iMax[1:(seasLag-1)]+1))
                    }
                }
            }
            # Start the loop for differences
            for(d in 1:nrow(iOrders)){
                m <- m + 1;
                if(!silent){
                    cat(paste0(rep("\b",nchar(round(m/nModelsARIMA,2)*100)+1),collapse=""));
                    cat(paste0(round((m)/nModelsARIMA,2)*100,"%"));
                }
                # If differences are zero, skip this step
                if(!all(iOrders[d,]==0)){
                    # Run the model for differences
                    testModel <- adam(y=yInSample, model=model, lags=lags,
                                      orders=list(ar=0,i=iOrders[d,],ma=0),
                                      distribution=distribution,
                                      h=h, holdout=FALSE,
                                      persistence=persistence, phi=phi, initial=initial,
                                      occurrence=occurrence, ic=ic, bounds=bounds,
                                      xreg=xreg, xregDo=xregDo, silent=TRUE, ...);
                }
                # Extract Information criteria
                ICValue <- ICFunction(testModel);
                if(silentDebug){
                    cat("I: "); cat(iOrders[d,]); cat(", "); cat(ICValue); cat("\n");
                }
                if(ICValue < bestICI){
                    bestICI <- ICValue;
                    dataMA <- dataI <- residuals(testModel);
                    if(ICValue < bestIC){
                        iBest <- iOrders[d,];
                        bestIC <- ICValue;
                        maBest <- arBest <- rep(0,length(arTest));
                    }
                }
                else{
                    if(fast){
                        m <- m + sum(maMax*(1 + sum(arMax)));
                        next;
                    }
                    else{
                        dataMA <- dataI <- residuals(testModel);
                    }
                }

                ##### Loop for MA #####
                if(any(maMax!=0)){
                    bestICMA <- bestICI;
                    maBestLocal <- maTest <- rep(0,length(maTest));
                    for(seasSelectMA in 1:length(lags)){
                        if(maMax[seasSelectMA]!=0){
                            for(maSelect in 1:maMax[seasSelectMA]){
                                m <- m + 1;
                                if(!silent){
                                    cat(paste0(rep("\b",nchar(round(m/nModelsARIMA,2)*100)+1),collapse=""));
                                    cat(paste0(round((m)/nModelsARIMA,2)*100,"%"));
                                }
                                maTest[seasSelectMA] <- maMax[seasSelectMA] - maSelect + 1;

                                # Run the model for MA
                                testModel <- adam(y=dataI, model="NNN", lags=lags,
                                                  orders=list(ar=0,i=0,ma=maTest),
                                                  distribution=distribution,
                                                  h=h, holdout=FALSE,
                                                  persistence=NULL, phi=NULL, initial=initial,
                                                  occurrence="none", ic=ic, bounds=bounds,
                                                  xreg=NULL, xregDo="use", silent=TRUE, ...);
                                # Exclude the variance from the number of parameters
                                nParamMA <- nparam(testModel)-1;
                                nParamNew <- nParamOriginal + nParamMA;
                                ICValue <- icCorrector(ICFunction(testModel), ic, nParamMA, obsNonzero, nParamNew);
                                if(silentDebug){
                                    cat("MA: "); cat(maTest); cat(", "); cat(ICValue); cat("\n");
                                }
                                if(ICValue < bestICMA){
                                    bestICMA <- ICValue;
                                    maBestLocal <- maTest;
                                    if(ICValue < bestIC){
                                        bestIC <- bestICMA;
                                        iBest <- iOrders[d,];
                                        maBest <- maTest;
                                        arBest <- rep(0,length(arTest));
                                    }
                                    dataMA <- residuals(testModel);
                                }
                                else{
                                    if(fast){
                                        m <- m + maTest[seasSelectMA] * (1 + sum(arMax)) - 1;
                                        maTest <- maBestLocal;
                                        break;
                                    }
                                    else{
                                        maTest <- maBestLocal;
                                        dataMA <- residuals(testModel);
                                    }
                                }

                                ##### Loop for AR #####
                                if(any(arMax!=0)){
                                    bestICAR <- bestICMA;
                                    arBestLocal <- arTest <- rep(0,length(arTest));
                                    for(seasSelectAR in 1:length(lags)){
                                        lagsTest[seasSelectAR] <- lags[seasSelectAR];
                                        if(arMax[seasSelectAR]!=0){
                                            for(arSelect in 1:arMax[seasSelectAR]){
                                                m <- m + 1;
                                                if(!silent){
                                                    cat(paste0(rep("\b",nchar(round(m/nModelsARIMA,2)*100)+1),collapse=""));
                                                    cat(paste0(round((m)/nModelsARIMA,2)*100,"%"));
                                                }
                                                arTest[seasSelectAR] <- arMax[seasSelectAR] - arSelect + 1;

                                                # Run the model for AR
                                                testModel <- adam(y=dataMA, model="NNN", lags=lags,
                                                                  orders=list(ar=arTest,i=0,ma=0),
                                                                  distribution=distribution,
                                                                  h=h, holdout=FALSE,
                                                                  persistence=NULL, phi=NULL, initial=initial,
                                                                  occurrence="none", ic=ic, bounds=bounds,
                                                                  xreg=NULL, xregDo="use", silent=TRUE, ...);
                                                # Exclude the variance from the number of parameters
                                                nParamAR <- nparam(testModel)-1;
                                                nParamNew <- nParamOriginal + nParamMA + nParamAR;
                                                ICValue <- icCorrector(ICFunction(testModel), ic, nParamAR, obsNonzero, nParamNew);
                                                if(silentDebug){
                                                    cat("AR: "); cat(arTest); cat(", "); cat(ICValue); cat("\n");
                                                }
                                                if(ICValue < bestICAR){
                                                    bestICAR <- ICValue;
                                                    arBestLocal <- arTest;
                                                    if(ICValue < bestIC){
                                                        bestIC <- ICValue;
                                                        iBest <- iOrders[d,];
                                                        arBest <- arTest;
                                                        maBest <- maTest;
                                                    }
                                                }
                                                else{
                                                    if(fast){
                                                        m <- m + arTest[seasSelectAR] - 1;
                                                        arTest <- arBestLocal;
                                                        break;
                                                    }
                                                    else{
                                                        arTest <- arBestLocal;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else{
                    ##### Loop for AR #####
                    if(any(arMax!=0)){
                        bestICAR <- bestICMA;
                        arBestLocal <- arTest <- rep(0,length(arTest));
                        for(seasSelectAR in 1:length(lags)){
                            lagsTest[seasSelectAR] <- lags[seasSelectAR];
                            if(arMax[seasSelectAR]!=0){
                                for(arSelect in 1:arMax[seasSelectAR]){
                                    m <- m + 1;
                                    if(!silent){
                                        cat(paste0(rep("\b",nchar(round(m/nModelsARIMA,2)*100)+1),collapse=""));
                                        cat(paste0(round((m)/nModelsARIMA,2)*100,"%"));
                                    }
                                    arTest[seasSelectAR] <- arMax[seasSelectAR] - arSelect + 1;
                                    nParamAR <- sum(arTest);
                                    nParamNew <- nParamOriginal + nParamAR;

                                    # Run the model for MA
                                    testModel <- adam(y=dataI, model="NNN", lags=lags,
                                                      orders=list(ar=arTest,i=0,ma=0),
                                                      distribution=distribution,
                                                      h=h, holdout=FALSE,
                                                      persistence=NULL, phi=NULL, initial=initial,
                                                      occurrence="none", ic=ic, bounds=bounds,
                                                      xreg=NULL, xregDo="use", silent=TRUE, ...);
                                    ICValue <- icCorrector(ICFunction(testModel), ic, nParamAR, obsNonzero, nParamNew);
                                    if(silentDebug){
                                        cat("AR: "); cat(arTest); cat(", "); cat(ICValue); cat("\n");
                                    }
                                    if(ICValue < bestICAR){
                                        bestICAR <- ICValue;
                                        arBestLocal <- arTest;
                                        if(ICValue < bestIC){
                                            bestIC <- ICValue;
                                            iBest <- iOrders[d,];
                                            arBest <- arTest;
                                            maBest <- maTest;
                                        }
                                    }
                                    else{
                                        if(fast){
                                            m <- m + arTest[seasSelectAR] - 1;
                                            arTest <- arBestLocal;
                                            break;
                                        }
                                        else{
                                            arTest <- arBestLocal;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(!silent && fast){
                cat(paste0(rep("\b",nchar(round(m/nModels,2)*100)+1),collapse=""));
                cat(paste0(" ",100,"%"));
            }

            #### Reestimate the best model in order to get rid of bias ####
            # Run the model for MA
            bestModel <- adam(y=y, model=modelOriginal, lags=lags,
                              orders=list(ar=(arBest),i=(iBest),ma=(maBest)),
                              distribution=distribution,
                              h=h, holdout=holdout,
                              persistence=persistenceOriginal, phi=phiOriginal, initial=initial,
                              occurrence=occurrenceOriginal, ic=ic, bounds=bounds,
                              xreg=xregOriginal, xregDo="use", silent=TRUE, ...);

            if(!silent){
                cat(". The best ARIMA is selected.\n");
            }
            return(bestModel);
        }
    }

    #### A simple loop, no ARIMA orders selection ####
    if(!arimaModelSelect){
        selectedModels <- adamReturner(y, model, lags, orders,
                                       distribution, h, holdout,
                                       persistence, phi, initial, arma,
                                       occurrence, ic, bounds,
                                       xreg, xregDo, parallel,
                                       arimaModelSelect, arMax, iMax, maMax, ...);
    }
    else{
        #### If there is ETS(X), do ARIMA selection on residuals ####
        # Extract residuals from adams for each distribution, fit best ARIMA for each, refit the models.
        if(etsModel || xregModel){
            selectedModels <- adamReturner(y, model, lags, orders,
                                           distribution, h, holdout,
                                           persistence, phi, initial, arma,
                                           occurrence, ic, bounds,
                                           xreg, xregDo, parallel,
                                           arimaModelSelect, arMax, iMax, maMax, ...);
        }
        #### Otherwise, do the stuff directly ####
        # Do ARIMA selection for each distribution in parallel.
        else{
            if(!parallel){
                # Prepare the list of models
                selectedModels <- vector("list",length(distribution));
                for(i in 1:length(distribution)){
                    if(!silent){
                        cat(paste0(distribution[i],": "));
                    }
                    selectedModels[[i]] <- arimaSelector(y=y, model=model,
                                                         lags=lags, arMax=arMax, iMax=iMax, maMax=maMax,
                                                         distribution=distribution[i], h=h, holdout=holdout,
                                                         persistence=persistence, phi=phi, initial=initial,
                                                         occurrence=occurrence, ic=ic, bounds=bounds, fast=fast,
                                                         silent=silent, xreg=xreg, xregDo=xregDo, testModelETS=NULL, ...);
                }
            }
            else{
                selectedModels <- foreach::`%dopar%`(foreach::foreach(i=1:length(distribution)),{
                    testModel <- arimaSelector(y=y, model=model,
                                               lags=lags, arMax=arMax, iMax=iMax, maMax=maMax,
                                               distribution=distribution[i], h=h, holdout=holdout,
                                               persistence=persistence, phi=phi, initial=initial,
                                               occurrence=occurrence, ic=ic, bounds=bounds, fast=fast,
                                               silent=TRUE, xreg=xreg, xregDo=xregDo, testModelETS=NULL, ...);
                    return(testModel);
                })
            }
        }
    }

    if(!silent){
        cat("Done!\n");
    }

    if(modelDo=="select"){
        ICValues <- sapply(selectedModels, ICFunction);
    }
    else{
        ICValues <- vector("numeric",length(distribution));
        for(i in 1:length(distribution)){
            ICValues[i] <- selectedModels[[i]]$ICs %*% selectedModels[[i]]$ICw;
        }
    }

    if(!silent){
        plot(selectedModels[[which.min(ICValues)]],7);
    }
    selectedModels[[which.min(ICValues)]]$timeElapsed <- Sys.time()-startTime;
    # names(ICValues) <- sapply(selectedModels, modelType);
    # selectedModels[[which.min(ICValues)]]$ICValues <- ICValues;

    return(selectedModels[[which.min(ICValues)]]);
}
