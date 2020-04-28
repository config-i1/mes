#' @rdname adam
#' @export
auto.adam <- function(y, model="ZXZ", lags=c(frequency(y)),
                      distribution=c("dnorm","dlogis","dlaplace","ds",
                                     "dlnorm","dllaplace","dls","dinvgauss"),
                      h=0, holdout=FALSE,
                      persistence=NULL, phi=NULL, initial=c("optimal","backcasting"),
                      occurrence=c("none","auto","fixed","general","odds-ratio","inverse-odds-ratio","direct"),
                      ic=c("AICc","AIC","BIC","BICc"), bounds=c("usual","admissible","none"),
                      xreg=NULL, xregDo=c("use","select"), xregInitial=NULL, xregPersistence=0,
                      silent=TRUE, ...){
    # Copyright (C) 2020 - Inf  Ivan Svetunkov

    # Start measuring the time of calculations
    startTime <- Sys.time();

    if(any(unlist(strsplit(model,""))=="C")){
        modelDo <- "combine";
    }
    else{
        modelDo <- "select";
    }

    selectedModels <- vector("list",length(distribution));
    if(!silent){
        cat("Evaluating models with different distributions... ");
    }
    for(i in 1:length(distribution)){
        if(!silent){
            cat(paste0(distribution[i],", "));
        }
        selectedModels[[i]] <- adam(y=y, model=model, lags=lags,
                                   distribution=distribution[i], loss="likelihood",
                                   h=h, holdout=holdout,
                                   persistence=persistence, phi=phi, initial=initial,
                                   occurrence=occurrence, ic=ic, bounds=bounds,
                                   xreg=xreg, xregDo=xregDo, xregInitial=xregInitial,
                                   xregPersistence=xregPersistence, silent=TRUE, ...);
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
