#' @param parallel If TRUE, the estimation of ADAM models is done in parallel (used in \code{auto.adam} only).
#' If the number is provided (e.g. \code{parallel=41}), then the specified number of cores is set up.
#' WARNING! Packages \code{foreach} and either \code{doMC} (Linux and Mac only)
#' or \code{doParallel} are needed in order to run the function in parallel.
#' @rdname adam
#' @export
auto.adam <- function(y, model="ZXZ", lags=c(frequency(y)), orders=list(ar=c(0),i=c(0),ma=c(0),select=FALSE),
                      distribution=c("dnorm","dlogis","dlaplace","ds",
                                     "dlnorm","dllaplace","dls","dinvgauss"),
                      h=0, holdout=FALSE,
                      persistence=NULL, phi=NULL, initial=c("optimal","backcasting"), arma=NULL,
                      occurrence=c("none","auto","fixed","general","odds-ratio","inverse-odds-ratio","direct"),
                      ic=c("AICc","AIC","BIC","BICc"), bounds=c("usual","admissible","none"),
                      xreg=NULL, xregDo=c("use","select","adapt"), silent=TRUE, parallel=FALSE, ...){
    # Copyright (C) 2020 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    if(any(unlist(strsplit(model,""))=="C")){
        modelDo <- "combine";
    }
    else{
        modelDo <- "select";
    }
    nModels <- length(distribution);

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

    if(!silent && !parallel){
        cat("Evaluating models with different distributions... ");
    }
    else if(!silent & parallel){
        cat(paste0("Working..."));
    }

    if(!parallel){
        # Prepare the list of models
        selectedModels <- vector("list",length(distribution));
        for(i in 1:length(distribution)){
            if(!silent){
                cat(paste0(distribution[i],", "));
            }
            selectedModels[[i]] <- adam(y=y, model=model, lags=lags, orders=orders,
                                        distribution=distribution[i],
                                        h=h, holdout=holdout,
                                        persistence=persistence, phi=phi, initial=initial, arma=arma,
                                        occurrence=occurrence, ic=ic, bounds=bounds,
                                        xreg=xreg, xregDo=xregDo, silent=TRUE, ...);
        }
    }
    else{
        selectedModels <- foreach(i=1:length(distribution)) %dopar% {
            return(adam(y=y, model=model, lags=lags, orders=orders,
                        distribution=distribution[i],
                        h=h, holdout=holdout,
                        persistence=persistence, phi=phi, initial=initial, arma=arma,
                        occurrence=occurrence, ic=ic, bounds=bounds,
                        xreg=xreg, xregDo=xregDo, silent=TRUE, ...));
        }
    }

    if(!silent){
        cat("Done!\n");
    }

    ic <- match.arg(ic,c("AICc","AIC","BIC","BICc"));
    ICFunction <- switch(ic,
                         "AIC"=AIC,
                         "AICc"=AICc,
                         "BIC"=BIC,
                         "BICc"=BICc);

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
