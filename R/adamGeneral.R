parametersChecker <- function(y, model, lags, formulaProvided, orders,
                              persistence, phi, initial,
                              distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                                             "dlnorm","dllaplace","dls","dinvgauss"),
                              loss, h, holdout,occurrence,
                              ic=c("AICc","AIC","BIC","BICc"), bounds=c("traditional","admissible","none"),
                              xreg, xregDo, responseName,
                              silent, modelDo, ParentEnvironment,
                              ellipsis, fast=FALSE){

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
        if(modelDo!="use"){
            lags <- frequency(y$x);
        }
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

    # Substitute NAs with mean values.
    yNAValues <- is.na(y);
    if(any(yNAValues)){
        warning("Data contains NAs. The values will be ignored during the model construction.",call.=FALSE);
        y[yNAValues] <- na.interp(y)[yNAValues];
    }

    # Define obs, the number of observations of in-sample
    obsAll <- length(y) + (1 - holdout)*h;
    obsInSample <- length(y) - holdout*h;

    # If this is just a numeric variable, use ts class
    if(all(yClasses=="integer") || all(yClasses=="data.frame") || all(yClasses=="matrix")){
        if(any(class(yIndex) %in% c("POSIXct","Date"))){
            yClasses <- "zoo";
        }
        else{
            yClasses <- "ts";
        }
    }
    yFrequency <- frequency(y);
    yStart <- yIndex[1];
    yInSample <- matrix(y[1:obsInSample],ncol=1);
    if(holdout){
        yForecastStart <- yIndex[obsInSample+1];
        yHoldout <- y[-c(1:obsInSample)];
        yForecastIndex <- yIndex[-c(1:obsInSample)];
        yInSampleIndex <- yIndex[c(1:obsInSample)];
    }
    else{
        yForecastStart <- yIndex[obsInSample]+as.numeric(diff(tail(yIndex,2)));
        yInSampleIndex <- yIndex;
        yForecastIndex <- yIndex[obsInSample]+as.numeric(diff(tail(yIndex,2)))*c(1:max(h,1));
        yHoldout <- NULL;
    }

    if(!is.numeric(yInSample)){
        stop("The provided data is not numeric! Can't construct any model!", call.=FALSE);
    }

    # Number of parameters to estimate / provided
    parametersNumber <- matrix(0,2,4,
                               dimnames=list(c("Estimated","Provided"),
                                             c("nParamInternal","nParamXreg","nParamOccurrence","nParamAll")));

    #### Check what is used for the model ####
    if(!is.character(model)){
        stop(paste0("Something strange is provided instead of character object in model: ",
                    paste0(model,collapse=",")),call.=FALSE);
    }

    # Predefine models pool for a model selection
    modelsPool <- NULL;
    if(!fast){
        # Deal with the list of models. Check what has been provided. Stop if there is a mistake.
        if(length(model)>1){
            if(any(nchar(model)>4)){
                stop(paste0("You have defined strange model(s) in the pool: ",
                            paste0(model[nchar(model)>4],collapse=",")),call.=FALSE);
            }
            else if(any(substr(model,1,1)!="A" & substr(model,1,1)!="M" & substr(model,1,1)!="C")){
                stop(paste0("You have defined strange model(s) in the pool: ",
                            paste0(model[substr(model,1,1)!="A" & substr(model,1,1)!="M"],collapse=",")),
                     call.=FALSE);
            }
            else if(any(substr(model,2,2)!="N" & substr(model,2,2)!="A" &
                        substr(model,2,2)!="M" & substr(model,2,2)!="C")){
                stop(paste0("You have defined strange model(s) in the pool: ",
                            paste0(model[substr(model,2,2)!="N" & substr(model,2,2)!="A" &
                                             substr(model,2,2)!="M"],collapse=",")),call.=FALSE);
            }
            else if(any(substr(model,3,3)!="N" & substr(model,3,3)!="A" &
                        substr(model,3,3)!="M" & substr(model,3,3)!="d" & substr(model,3,3)!="C")){
                stop(paste0("You have defined strange model(s) in the pool: ",
                            paste0(model[substr(model,3,3)!="N" & substr(model,3,3)!="A" &
                                             substr(model,3,3)!="M" & substr(model,3,3)!="d"],collapse=",")),
                     call.=FALSE);
            }
            else if(any(nchar(model)==4 & substr(model,4,4)!="N" & substr(model,4,4)!="A" &
                        substr(model,4,4)!="M" & substr(model,4,4)!="C")){
                stop(paste0("You have defined strange model(s) in the pool: ",
                            paste0(model[nchar(model)==4 & substr(model,4,4)!="N" &
                                             substr(model,4,4)!="A" & substr(model,4,4)!="M"],collapse=",")),
                     call.=FALSE);
            }
            else{
                modelsPoolCombiner <- (substr(model,1,1)=="C" | substr(model,2,2)=="C" |
                                           substr(model,3,3)=="C" | substr(model,4,4)=="C");
                modelsPool <- model[!modelsPoolCombiner];
                modelsPool <- unique(modelsPool);
                if(any(modelsPoolCombiner)){
                    if(any(substr(model,nchar(model),nchar(model))!="N")){
                        model <- "CCC";
                    }
                    else{
                        model <- "CCN";
                    }
                }
                else{
                    model <- c("Z","Z","Z");
                    if(all(substr(modelsPool,nchar(modelsPool),nchar(modelsPool))=="N")){
                        model[3] <- "N";
                    }
                    if(all(substr(modelsPool,2,2)=="N")){
                        model[2] <- "N";
                    }
                    model <- paste0(model,collapse="");
                }
            }
        }
    }

    # If chosen model is "AAdN" or anything like that, we are taking the appropriate values
    if(nchar(model)==4){
        Etype <- substr(model,1,1);
        Ttype <- substr(model,2,2);
        Stype <- substr(model,4,4);
        damped <- TRUE;
        if(substr(model,3,3)!="d"){
            message(paste0("You have defined a strange model: ",model));
            model <- paste0(Etype,Ttype,"d",Stype);
        }
    }
    else if(nchar(model)==3){
        Etype <- substr(model,1,1);
        Ttype <- substr(model,2,2);
        Stype <- substr(model,3,3);
        if(any(Ttype==c("Z","X","Y"))){
            damped <- TRUE;
        }
        else{
            damped <- FALSE;
        }
    }
    else{
        message(paste0("You have defined a strange model: ",model));
        message("Switching to 'ZZZ'");
        model <- "ZZZ";

        Etype <- "Z";
        Ttype <- "Z";
        Stype <- "Z";
        damped <- TRUE;
    }

    # Define if we want to select or combine models... or do none of the above.
    if(is.null(modelsPool)){
        if(any(unlist(strsplit(model,""))=="C")){
            modelDo <- "combine";
            if(Etype=="C"){
                Etype <- "Z";
            }
            if(Ttype=="C"){
                Ttype <- "Z";
            }
            if(Stype=="C"){
                Stype <- "Z";
            }
        }
        else if(any(unlist(strsplit(model,""))=="Z") ||
                any(unlist(strsplit(model,""))=="X") ||
                any(unlist(strsplit(model,""))=="Y") ||
                any(unlist(strsplit(model,""))=="F") ||
                any(unlist(strsplit(model,""))=="P")){
            modelDo <- "select";

            # The full test, sidestepping branch and bound
            if(any(unlist(strsplit(model,""))=="F")){
                modelsPool <- c("ANN","AAN","AAdN","AMN","AMdN",
                                "ANA","AAA","AAdA","AMA","AMdA",
                                "ANM","AAM","AAdM","AMM","AMdM",
                                "MNN","MAN","MAdN","MMN","MMdN",
                                "MNA","MAA","MAdA","MMA","MMdA",
                                "MNM","MAM","MAdM","MMM","MMdM");
                Etype[] <- Ttype[] <- Stype[] <- "Z";
                model <- "FFF";
            }

            # The test for pure models only
            if(any(unlist(strsplit(model,""))=="P")){
                modelsPool <- c("ANN","AAN","AAdN","ANA","AAA","AAdA",
                                "MNN","MMN","MMdN","MNM","MMM","MMdM");
                Etype[] <- Ttype[] <- Stype[] <- "Z";
                model <- "PPP";
            }
        }
        else{
            modelDo <- "estimate";
        }

        if(Etype=="X"){
            Etype <- "A";
        }
        else if(Etype=="Y"){
            Etype <- "M";
        }
    }
    else{
        if(any(unlist(strsplit(model,""))=="C")){
            modelDo <- "combine";
        }
        else{
            modelDo <- "select";
        }
    }

    modelIsTrendy <- (Ttype!="N");

    #### Check the components of model ####
    componentsNamesETS <- "level";
    componentsNumberETS <- 1;
    ### Check error type
    if(all(Etype!=c("Z","X","Y","A","M","C"))){
        warning(paste0("Wrong error type: ",Etype,". Should be 'Z', 'X', 'Y', 'A' or 'M'. ",
                       "Changing to 'Z'"),call.=FALSE);
        Etype <- "Z";
        modelDo <- "select";
    }

    ### Check trend type
    if(all(Ttype!=c("Z","X","Y","N","A","M","C"))){
        warning(paste0("Wrong trend type: ",Ttype,". Should be 'Z', 'X', 'Y', 'N', 'A' or 'M'. ",
                       "Changing to 'Z'"),call.=FALSE);
        Ttype <- "Z";
        modelDo <- "select";
    }
    if(modelIsTrendy){
        componentsNamesETS <- c(componentsNamesETS,"trend");
        componentsNumberETS[] <- componentsNumberETS+1;
    }

    #### Check the lags vector ####
    if(any(c(lags)<0)){
        stop("Right! Why don't you try complex lags then, mister smart guy?",call.=FALSE);
    }

    # If there are zero lags, drop them
    if(any(lags==0)){
        lags <- lags[lags!=0];
    }

    # Form the lags based on the provided stuff. Get rid of ones and leave unique seasonals
    # Add one for the level
    lags <- c(1,unique(lags[lags>1]));

    #### ARIMA term ####
    # This should be available for pure models only
    if(is.list(orders)){
        arOrders <- orders$ar;
        iOrders <- orders$i;
        maOrders <- orders$ma;
    }
    else if(is.vector(orders)){
        arOrders <- orders[1];
        iOrders <- orders[2];
        maOrders <- orders[3];
    }

    # If there is arima, prepare orders
    if(sum(c(arOrders,iOrders,maOrders))>0){
        arimaModel <- TRUE;

        # See if AR is needed
        arRequired <- FALSE;
        if(sum(arOrders)>0){
            arRequired[] <- TRUE;
        }

        # See if I is needed
        iRequired <- FALSE;
        if(sum(iOrders)>0){
            iRequired[] <- TRUE;
        }

        # See if I is needed
        maRequired <- FALSE;
        if(sum(maOrders)>0){
            maRequired[] <- TRUE;
        }

        # Define maxOrder and make all the values look similar (for the polynomials)
        maxOrder <- max(length(arOrders),length(iOrders),length(maOrders),length(lags));
        if(length(arOrders)!=maxOrder){
            arOrders <- c(arOrders,rep(0,maxOrder-length(arOrders)));
        }
        if(length(iOrders)!=maxOrder){
            iOrders <- c(iOrders,rep(0,maxOrder-length(iOrders)));
        }
        if(length(maOrders)!=maxOrder){
            maOrders <- c(maOrders,rep(0,maxOrder-length(maOrders)));
        }

        # Define the non-zero values. This is done via the calculation of orders of polynomials
        ariValues <- list(NA);
        maValues <- list(NA);
        for(i in 1:length(lags)){
            ariValues[[i]] <- c(0,min(1,arOrders[i]):arOrders[i])
            if(iOrders[i]!=0){
                ariValues[[i]] <- c(ariValues[[i]],1:iOrders[i]+arOrders[i]);
            }
            ariValues[[i]] <- unique(ariValues[[i]] * lags[i]);
            maValues[[i]] <- unique(c(0,min(1,maOrders[i]):maOrders[i]) * lags[i]);
        }

        # Produce ARI polynomials
        ariLengths <- unlist(lapply(ariValues,length));
        ariPolynomial <- array(0,ariLengths);
        for(i in 1:length(ariValues)){
            if(i==1){
                ariPolynomial <- ariPolynomial + array(ariValues[[i]], ariLengths);
            }
            else{
                ariPolynomial <- ariPolynomial + array(rep(ariValues[[i]],each=prod(ariLengths[1:(i-1)])),
                                                       ariLengths);
            }
        }

        # Produce MA polynomials
        maLengths <- unlist(lapply(maValues,length));
        maPolynomial <- array(0,maLengths);
        for(i in 1:length(maValues)){
            if(i==1){
                maPolynomial <- maPolynomial + array(maValues[[i]], maLengths);
            }
            else{
                maPolynomial <- maPolynomial + array(rep(maValues[[i]],each=prod(maLengths[1:(i-1)])),
                                                     maLengths);
            }
        }

        # What are the non-zero ARI and MA polynomials?
        ### What are their positions in transition matrix?
        nonZeroARI <- unique(matrix(c(ariPolynomial)[-1],ncol=1));
        nonZeroMA <- unique(matrix(c(maPolynomial)[-1],ncol=1));
        # Lags for the ARIMA components
        lagsModelARIMA <- matrix(sort(unique(c(nonZeroARI,nonZeroMA))),ncol=1);
        nonZeroARI <- cbind(nonZeroARI+1,which(lagsModelARIMA %in% nonZeroARI));
        nonZeroMA <- cbind(nonZeroMA+1,which(lagsModelARIMA %in% nonZeroMA));

        # Number of components
        componentsNumberARIMA <- length(lagsModelARIMA);
        # Their names
        componentsNamesARIMA <- paste0("ARIMAState",c(1:componentsNumberARIMA));
        # Number of initials needed. This is based on the longest one. The others are just its transformations
        initialArimaNumber <- max(lagsModelARIMA);

        if(obsInSample < initialArimaNumber){
            warning(paste0("In-sample size is ",obsInSample,", while number of ARIMA components is ",componentsNumberARIMA,
                           ". Cannot fit the model."),call.=FALSE)
            stop("Not enough observations for such a complicated model.",call.=FALSE);
        }
    }
    else{
        arOrders <- NULL;
        iOrders <- NULL;
        maOrders <- NULL;
        arimaModel <- FALSE;
        arRequired <- arEstimate <- FALSE;
        iRequired <- FALSE;
        maRequired <- maEstimate <- FALSE;
        lagsModelARIMA <- initialArimaNumber <- 0;
        componentsNumberARIMA <- 0;
        componentsNamesARIMA <- NULL;
        nonZeroARI <- NULL;
        nonZeroMA <- NULL;
    }

    modelIsSeasonal <- Stype!="N";

    # Lags of the model used inside the functions
    lagsModel <- matrix(lags,ncol=1);

    # If we have a trend add one more lag
    if(modelIsTrendy){
        lagsModel <- rbind(1,lagsModel);
    }
    # If we don't have seasonality, remove seasonal lag
    if(!modelIsSeasonal & any(lagsModel>1)){
        lagsModel <- lagsModel[lagsModel==1,,drop=FALSE];
    }

    # Lags of the model
    lagsModelSeasonal <- lagsModel[lagsModel>1];
    lagsModelMax <- max(lagsModel);
    lagsLength <- length(lagsModel);

    #### Check the seasonal model vs lags ####
    if(all(Stype!=c("Z","X","Y","N","A","M","C"))){
        warning(paste0("Wrong seasonality type: ",Stype,". Should be 'Z', 'X', 'Y', 'C', 'N', 'A' or 'M'. ",
                       "Setting to 'Z'."),call.=FALSE);
        if(lagsModelMax==1){
            Stype <- "N";
            modelIsSeasonal <- FALSE;
        }
        else{
            Stype <- "Z";
            modelDo <- "select";
        }
    }
    if(all(modelIsSeasonal,lagsModelMax==1)){
        if(all(Stype!=c("Z","X","Y"))){
            warning(paste0("Cannot build the seasonal model on data with the unity lags.\n",
                           "Switching to non-seasonal model: ETS(",substr(model,1,nchar(model)-1),"N)"));
        }
        Stype <- "N";
        modelIsSeasonal <- FALSE;
        substr(model,nchar(model),nchar(model)) <- "N";
    }

    # Check the pool of models to combine if it was decided that the data is not seasonal
    if(!modelIsSeasonal && !is.null(modelsPool)){
        modelsPool <- modelsPool[substr(modelsPool,nchar(modelsPool),nchar(modelsPool))=="N"];
    }

    # Check the type of seasonal
    if(Stype!="N"){
        componentsNamesETS <- c(componentsNamesETS,"seasonal");
        componentsNumberETS[] <- componentsNumberETS+1;
        componentsNumberETSSeasonal <- 1;
    }
    else{
        componentsNumberETSSeasonal <- 0;
    }

    # Check, whether the number of lags and the number of components are the same
    if(lagsLength>componentsNumberETS){
        if(Stype!="N"){
            componentsNamesETS <- c(componentsNamesETS[-length(componentsNamesETS)],paste0("seasonal",c(1:(lagsLength-componentsNumberETS-1))));
            componentsNumberETSSeasonal[] <- lagsLength-componentsNumberETS+1;
            componentsNumberETS[] <- lagsLength;
        }
        else{
            lagsModel <- lagsModel[1:componentsNumberETS,,drop=FALSE];
            lagsModelMax <- max(lagsModel);
            lagsLength <- length(lagsModel);
        }
    }
    else if(lagsLength<componentsNumberETS){
        stop("The number of components of the model is smaller than the number of provided lags", call.=FALSE);
    }

    if(!fast){
        #### Distribution selected ####
        distribution <- match.arg(distribution);
    }

    #### Loss function type ####
    if(is.function(loss)){
        lossFunction <- loss;
        loss <- "custom";
        multisteps <- FALSE;
    }
    else{
        loss <- match.arg(loss[1],c("likelihood","MSE","MAE","HAM","LASSO","RIDGE",
                                    "MSEh","TMSE","GTMSE","MSCE",
                                    "MAEh","TMAE","GTMAE","MACE",
                                    "HAMh","THAM","GTHAM","CHAM","GPL",
                                    "aMSEh","aTMSE","aGTMSE","aMSCE","aGPL"));

        if(any(loss==c("MSEh","TMSE","GTMSE","MSCE","MAEh","TMAE","GTMAE","MACE",
                       "HAMh","THAM","GTHAM","CHAM","GPL",
                       "aMSEh","aTMSE","aGTMSE","aMSCE","aGPL"))){
            if(!is.null(h) && h>0){
                multisteps <- TRUE;
            }
            else{
                stop("The horizon \"h\" needs to be specified and be positive in order for the multistep loss to work.",
                     call.=FALSE);
                multisteps <- FALSE;
            }
        }
        else{
            multisteps <- FALSE;
        }
        if(any(loss==c("LASSO","RIDGE"))){
            warning(paste0(loss," is not yet implemented properly. This is an experimental option. Use with care."),
                    call.=FALSE);
        }
        lossFunction <- NULL;
    }

    #### Explanatory variables: xregExist and xregDo ####
    xregDo <- match.arg(xregDo,c("use","select","adapt"));
    xregExist <- !is.null(xreg);

    #### Persistence provided ####
    # Vectors for persistence of different components
    persistenceLevel <- NULL;
    persistenceTrend <- NULL;
    persistenceSeasonal <- NULL;
    persistenceXreg <- NULL;
    # InitialEstimate vectors, defining what needs to be estimated
    persistenceEstimate <- persistenceLevelEstimate <- persistenceTrendEstimate <-
        persistenceXregEstimate <- TRUE;
    # persistence of seasonal is a vector, not a scalar, because we can have several lags
    persistenceSeasonalEstimate <- rep(TRUE,componentsNumberETSSeasonal);
    if(!is.null(persistence)){
        # If it is a list
        if(is.list(persistence)){
            # If this is a named list, then extract stuff using names
            if(!is.null(names(persistence))){
                if(!is.null(persistence$level)){
                    persistenceLevel <- persistence$level;
                }
                else if(!is.null(persistence$alpha)){
                    persistenceLevel <- persistence$alpha;
                }
                if(!is.null(persistence$trend)){
                    persistenceTrend <- persistence$trend;
                }
                else if(!is.null(persistence$beta)){
                    persistenceTrend <- persistence$beta;
                }
                if(!is.null(persistence$seasonal)){
                    persistenceSeasonal <- persistence$seasonal;
                }
                else if(!is.null(persistence$gamma)){
                    persistenceSeasonal <- persistence$gamma;
                }
                if(!is.null(persistence$xreg)){
                    persistenceXreg <- persistence$xreg;
                }
                else if(!is.null(persistence$delta)){
                    persistenceXreg <- persistence$delta;
                }
            }
            else{
                if(!is.null(persistence[[1]])){
                    persistenceLevel <- persistence[[1]];
                }
                if(!is.null(persistence[[2]])){
                    persistenceTrend <- persistence[[2]];
                }
                if(!is.null(persistence[[3]])){
                    persistenceSeasonal <- persistence[[3]];
                }
                if(!is.null(persistence[[4]])){
                    persistenceXreg <- persistence[[4]];
                }
            }
            # Define estimate variables
            if(!is.null(persistenceLevel)){
                persistenceLevelEstimate[] <- FALSE;
                parametersNumber[2,1] <- parametersNumber[2,1] + 1;
            }
            if(!is.null(persistenceTrend)){
                persistenceTrendEstimate[] <- FALSE;
                parametersNumber[2,1] <- parametersNumber[2,1] + 1;
            }
            if(!is.null(persistenceSeasonal)){
                if(is.list(persistenceSeasonal)){
                    persistenceSeasonalEstimate[] <- length(persistenceSeasonal)==length(lagsModelSeasonal);
                }
                else{
                    persistenceSeasonalEstimate[] <- FALSE;
                }
                parametersNumber[2,1] <- parametersNumber[2,1] + length(unlist(persistenceSeasonal));
            }
            if(!is.null(persistenceXreg)){
                persistenceXregEstimate[] <- FALSE;
                parametersNumber[2,1] <- parametersNumber[2,1] + length(persistenceXreg);
            }
        }
        else if(is.numeric(persistence)){
            if(modelDo!="estimate"){
                warning(paste0("Predefined persistence vector can only be used with ",
                               "preselected ETS model.\n",
                               "Changing to estimation of persistence vector values."),call.=FALSE);
                persistence <- NULL;
                persistenceEstimate <- TRUE;
            }
            else{
                # If it is smaller... We don't know the length of xreg yet at this stage
                if(length(persistence)<lagsLength){
                    warning(paste0("Length of persistence vector is wrong! ",
                                   "Changing to estimation of persistence vector values."),
                            call.=FALSE);
                    persistence <- NULL;
                    persistenceEstimate <- TRUE;
                }
                else{
                    j <- 1;
                    persistenceLevel <- as.vector(persistence)[1];
                    names(persistenceLevel) <- "alpha";
                    if(modelIsTrendy && length(persistence)>j){
                        j <- j+1;
                        persistenceTrend <- as.vector(persistence)[j];
                        names(persistenceTrend) <- "beta";
                    }
                    if(Stype!="N" && length(persistence)>j){
                        j <- j+1;
                        persistenceSeasonal <- as.vector(persistence)[j];
                        names(persistenceSeasonal) <- paste0("gamma",c(1:length(persistenceSeasonal)));
                    }
                    if(xregExist && length(persistence)>j){
                        persistenceXreg <- as.vector(persistence)[-c(1:j)];
                        names(persistenceXreg) <- paste0("delta",c(1:length(persistenceXreg)));
                    }

                    persistenceEstimate[] <- persistenceLevelEstimate[] <- persistenceTrendEstimate[] <-
                        persistenceXregEstimate[] <- persistenceSeasonalEstimate[] <- FALSE;
                    parametersNumber[2,1] <- parametersNumber[2,1] + length(persistence);
                    bounds <- "n";
                }
            }
        }
        else{
            warning(paste0("Persistence is not a numeric vector!\n",
                           "Changing to estimation of persistence vector values."),call.=FALSE);
            persistence <- NULL;
            persistenceEstimate <- TRUE;
        }
    }
    else{
        persistenceEstimate <- TRUE;
    }

   # Make sure that only important elements are estimated.
    if(Ttype=="N"){
        persistenceTrendEstimate[] <- FALSE;
        persistenceTrend <- NULL;
    }
    if(Stype=="N"){
        persistenceSeasonalEstimate[] <- FALSE;
        persistenceSeasonal <- NULL;
    }
    if(!xregExist){
        persistenceXregEstimate[] <- FALSE;
        persistenceXreg <- NULL;
    }

    # Redefine persitenceEstimate value
    persistenceEstimate[] <- any(c(persistenceLevelEstimate,persistenceTrendEstimate,
                                 persistenceSeasonalEstimate,persistenceXregEstimate));

    #### Phi ####
    if(!is.null(phi)){
        if(!is.numeric(phi) & (damped)){
            warning(paste0("Provided value of phi is meaningless. phi will be estimated."),
                    call.=FALSE);
            phi <- 0.95;
            phiEstimate <- TRUE;
        }
        else if(is.numeric(phi) & (phi<0 | phi>2)){
            warning(paste0("Damping parameter should lie in (0, 2) region. ",
                           "Changing to the estimation of phi."),call.=FALSE);
            phi[] <- 0.95;
            phiEstimate <- TRUE;
        }
        else{
            phiEstimate <- FALSE;
            if(damped){
                parametersNumber[2,1] <- parametersNumber[2,1] + 1;
            }
        }
    }
    else{
        if(damped){
            phiEstimate <- TRUE;
            phi <- 0.95;
        }
        else{
            phiEstimate <- FALSE;
            phi <- 1;
        }
    }


    #### Lags for ARIMA ####
    if(arimaModel){
        lagsModelAll <- rbind(lagsModel,lagsModelARIMA);
        lagsModelMax <- max(lagsModel);
    }
    else{
        lagsModelAll <- lagsModel;
    }


    #### Occurrence variable ####
    if(is.occurrence(occurrence)){
        oesModel <- occurrence;
        occurrence <- oesModel$occurrence;
        if(oesModel$occurrence=="provided"){
            occurrenceModelProvided <- FALSE;
        }
        else{
            occurrenceModelProvided <- TRUE;
        }
    }
    else{
        occurrenceModelProvided <- FALSE;
        oesModel <- NULL;
    }
    pFitted <- matrix(1, obsInSample, 1);
    pForecast <- rep(NA,h);

    if(is.numeric(occurrence)){
        # If it is data, then it should correspond to the in-sample.
        if(all(occurrence==1)){
            occurrence <- "none";
        }
        else{
            if(any(occurrence<0,occurrence>1)){
                warning(paste0("Parameter 'occurrence' should contain values between zero and one.\n",
                               "Converting to appropriate vector."),call.=FALSE);
                occurrence[] <- (occurrence!=0)*1;
            }

            # "provided", meaning that we have been provided the values of p
            pFitted[] <- occurrence[1:obsInSample];
            # Create forecasted values for occurrence
            if(h>0){
                if(length(occurrence)>obsInSample){
                    pForecast <- occurrence[-c(1:obsInSample)];
                }
                else{
                    pForecast <- rep(tail(occurrence,1),h);
                }
                if(length(pForecast)>h){
                    pForecast <- pForecast[1:h];
                }
                else if(length(pForecast)<h){
                    pForecast <- c(pForecast,rep(tail(pForecast,1),h-length(pForecast)));
                }
            }
            else{
                pForecast <- NA;
            }
            occurrence <- "provided";
            oesModel <- list(fitted=pFitted,forecast=pForecast,occurrence="provided");
        }
    }

    occurrence <- match.arg(occurrence[1],c("none","auto","fixed","general","odds-ratio",
                                            "inverse-odds-ratio","direct","provided"));

    otLogical <- yInSample!=0;

    # If the data is not occurrence, let's assume that the parameter was switched unintentionally.
    if(all(otLogical) & all(occurrence!=c("none","provided"))){
        occurrence <- "none";
        occurrenceModelProvided <- FALSE;
    }

    # If there were NAs and the occurrence was not specified, do something with it
    # In all the other cases, NAs will be sorted out by the model
    if(any(yNAValues) && (occurrence=="none")){
        otLogical <- (!yNAValues)[1:obsInSample];
        occurrence[] <- "provided";
        pFitted[] <- otLogical*1;
        pForecast[] <- 1;
        occurrenceModel <- FALSE;
        oesModel <- structure(list(y=matrix((otLogical)*1,ncol=1),fitted=pFitted,forecast=pForecast,
                                   occurrence="provided"),class="occurrence");
    }
    else{
        if(occurrence=="none"){
            occurrenceModel <- FALSE;
            otLogical <- rep(TRUE,obsInSample);
        }
        else if(occurrence=="provided"){
            occurrenceModel <- FALSE;
            oesModel$y <- matrix(otLogical*1,ncol=1);
        }
        else{
            occurrenceModel <- TRUE;
        }

        # Include NAs in the zeroes, so that we avoid weird cases
        if(any(yNAValues)){
            otLogical <- !(!otLogical | yNAValues[1:obsInSample]);
        }
    }

    if(any(yClasses=="ts")){
        ot <- ts(matrix(otLogical*1,ncol=1), start=yStart, frequency=yFrequency);
    }
    else{
        ot <- ts(matrix(otLogical*1,ncol=1), start=c(0,0), frequency=lagsModelMax);
    }
    obsNonzero <- sum(ot);
    obsZero <- obsInSample - obsNonzero;


    #### Initial values ####
    # Vectors for initials of different components
    initialLevel <- NULL;
    initialTrend <- NULL;
    initialSeasonal <- NULL;
    initialArima <- NULL;
    initialXreg <- NULL;
    # InitialEstimate vectors, defining what needs to be estimated
    # NOTE: that initial==c("optimal","backcasting") meanst initialEstimate==TRUE!
    initialEstimate <- initialLevelEstimate <- initialTrendEstimate <-
        initialArimaEstimate <- initialXregEstimate <- TRUE;
    # initials of seasonal is a vector, not a scalar, because we can have several lags
    initialSeasonalEstimate <- rep(TRUE,componentsNumberETSSeasonal);

    # This is an initialisation of the variable
    initialType <- "optimal"
    # initial type can be: "o" - optimal, "b" - backcasting, "p" - provided.
    if(any(is.character(initial))){
        initialType[] <- match.arg(initial, c("optimal","backcasting"));
    }
    else if(is.null(initial)){
        if(!silent){
            message("Initial value is not selected. Switching to optimal.");
        }
        initialType[] <- "optimal";
    }
    else if(!is.null(initial)){
        if(modelDo!="estimate"){
            warning(paste0("Predefined initials vector can only be used with preselected ETS model.\n",
                           "Changing to estimation of initials."),call.=FALSE);
            initialType[] <- "optimal";
            initialEstimate[] <- initialLevelEstimate[] <- initialTrendEstimate[] <-
                initialSeasonalEstimate[] <- initialArimaEstimate[] <- initialXregEstimate[] <- TRUE;
        }
        else{
            # If the list is provided, then check what this is.
            # This should be: level, trend, seasonal[[1]], seasonal[[2]], ..., ARIMA, xreg
            if(is.list(initial)){
                # If this is a named list, then extract stuff using names
                if(!is.null(names(initial))){
                    if(!is.null(initial$level)){
                        initialLevel <- initial$level;
                    }
                    if(!is.null(initial$trend)){
                        initialTrend <- initial$trend;
                    }
                    if(!is.null(initial$seasonal)){
                        initialSeasonal <- initial$seasonal;
                    }
                    if(!is.null(initial$ARIMA)){
                        initialArima <- initial$ARIMA;
                    }
                    if(!is.null(initial$xreg)){
                        initialXreg <- initial$xreg;
                    }
                }
                else{
                    if(!is.null(initial[[1]])){
                        initialLevel <- initial[[1]];
                    }
                    if(!is.null(initial[[2]])){
                        initialTrend <- initial[[2]];
                    }
                    if(!is.null(initial[[3]])){
                        initialSeasonal <- initial[[3]];
                    }
                    if(!is.null(initial[[4]])){
                        initialArima <- initial[[4]];
                    }
                    if(!is.null(initial[[5]])){
                        initialXreg <- initial[[5]];
                    }
                }
            }
            else{
                if(!is.numeric(initial)){
                    warning(paste0("Initial vector is not numeric!\n",
                                   "Values of initial vector will be estimated."),call.=FALSE);
                    initialType[] <- "optimal";
                }
                else{
                    # If this is a vector, then it should contain values in the order:
                    # level, trend, seasonal1, seasonal2, ..., ARIMA, xreg
                    # if(length(initial)<(sum(lagsModelAll))){
                    #     warning(paste0("The vector of initials contains only values for several components. ",
                    #                    "We will use what we can."),call.=FALSE);
                    # }
                    # else{
                        j <- 1;
                        initialLevel <- initial[1];
                        initialLevelEstimate[] <- FALSE;
                        if(modelIsTrendy){
                            j <- 2;
                            # If there is something in the vector, use it
                            if(all(!is.na(initial[j]))){
                                initialTrend <- initial[j];
                                initialTrendEstimate[] <- FALSE;
                            }
                        }
                        if(Stype!="N"){
                            # If there is something in the vector, use it
                            if(length(initial[-c(1:j)])>0){
                                initialSeasonal <- vector("list",componentsNumberETSSeasonal);
                                m <- 0;
                                for(i in 1:componentsNumberETSSeasonal){
                                    if(all(!is.na(initial[j+m+1:lagsModelSeasonal[i]]))){
                                        initialSeasonal[[i]] <- initial[j+m+1:lagsModelSeasonal[i]];
                                        m <- m + lagsModelSeasonal[i];
                                    }
                                    else{
                                        break;
                                    }
                                }
                                j <- j+m;
                                initialSeasonalEstimate[] <- FALSE;
                            }
                        }
                        if(arimaModel){
                            # If there is something else left, this must be ARIMA
                            if(all(!is.na(initial[j+c(1:initialArimaNumber)]))){
                                initialArima <- initial[j+c(1:initialArimaNumber)];
                                j <- j+max(lagsModelARIMA);
                                initialArimaEstimate[] <- FALSE;
                            }
                        }
                        if(xregExist){
                            # Something else? xreg for sure!
                            if(length(initial[-c(1:j)])>0){
                                initialXreg <- initial[-c(1:j)];
                                initialXregEstimate[] <- FALSE
                            }
                        }
                        parametersNumber[2,1] <- parametersNumber[2,1] + j;
                    }
                # }
            }
        }
    }

    #### Check the provided initials and define initialEstimate variables ####
    # Level
    if(!is.null(initialLevel)){
        if(length(initialLevel)>1){
            warning("Initial level contains more than one value! Using the first one.",
                    call.=FALSE);
            initialLevel <- initialLevel[1];
        }
        initialLevelEstimate[] <- FALSE;
        parametersNumber[2,1] <- parametersNumber[2,1] + 1;
    }
    # Trend
    if(!is.null(initialTrend)){
        if(length(initialTrend)>1){
            warning("Initial trend contains more than one value! Using the first one.",
                    call.=FALSE);
            initialTrend <- initialTrend[1];
        }
        initialTrendEstimate[] <- FALSE;
        parametersNumber[2,1] <- parametersNumber[2,1] + 1;
    }
    # Seasonal
    if(!is.null(initialSeasonal)){
        # The list means several seasonal lags
        if(is.list(initialSeasonal)){
            # Is the number of seasonal initials correct? If it is bigger, then remove redundant
            if(length(initialSeasonal)>componentsNumberETSSeasonal){
                warning("Initial seasonals contained more elements than needed! Removing redundant ones.",
                        call.=FALSE);
                initialSeasonal <- initialSeasonal[1:componentsNumberETSSeasonal];
            }
            # Is the number of initials in each season correct? Use the correct ones only
            if(any(!(sapply(initialSeasonal,length) %in% lagsModelSeasonal))){
                warning(paste0("Some of initial seasonals have a wrong length, ",
                               "not corresponding to the provided lags. We will estimate them."),
                        call.=FALSE);
                initialSeasonalToUse <- sapply(initialSeasonal,length) %in% lagsModelSeasonal;
                initialSeasonal <- initialSeasonal[initialSeasonalToUse];
            }
            initialSeasonalEstimate[] <- !(lagsModelSeasonal %in% sapply(initialSeasonal,length));
            # If there are some gaps in what to estimate, reform initialSeason to make sense in the future creator function
            if(!all(initialSeasonalEstimate) && !all(!initialSeasonalEstimate)){
                initialSeasonalCorrect <- vector("list",componentsNumberETSSeasonal);
                initialSeasonalCorrect[which(!initialSeasonalEstimate)] <- initialSeasonal;
                initialSeasonal <- initialSeasonalCorrect;
            }
        }
        # The vector implies only one seasonal
        else{
            if(all(length(initialSeasonal)!=lagsModelSeasonal)){
                warning(paste0("Wrong length of seasonal initial: ",length(initialSeasonal),
                               "Instead of ",lagsModelSeasonal,". Switching to estimation."),
                        call.=FALSE)
                initialSeasonalEstimate[] <- TRUE;
            }
            else{
                initialSeasonalEstimate[] <- FALSE;
            }
            # Create a list from the vector for consistency purposes
            initialSeasonal <- list(initialSeasonal);
        }
        parametersNumber[2,1] <- parametersNumber[2,1] + length(unlist(initialSeasonal));
    }
    # ARIMA
    if(!is.null(initialArima)){
        if(length(initialArima)!=initialArimaNumber){
            warning(paste0("The length of ARIMA initials is ",length(initialArima),
                           " instead of ",initialArimaNumber,". Estimating initials instead!"),
                    call.=FALSE);
            initialArimaEstimate[] <- TRUE;
        }
        else{
            initialArimaEstimate[] <- FALSE;
            parametersNumber[2,1] <- parametersNumber[2,1] + length(initialArima);
        }
    }
    # xreg
    if(!is.null(initialXreg)){
        initialXregEstimate[] <- FALSE;
        parametersNumber[2,1] <- parametersNumber[2,1] + length(initialXreg);
    }

    #### Check ARIMA parameters, if they are provided ####
    if(arimaModel){
        arEstimate <- maEstimate <- FALSE;
        if(any(arOrders>0)){
            arEstimate[] <- TRUE;
        }
        if(any(maOrders>0)){
            maEstimate[] <- TRUE;
        }
        # Check the provided parameters for AR and MA
    }
    else{}


    #### xreg preparation ####
    # Check the xregDo
    if(!xregExist){
        xregDo[] <- "use";
        formulaProvided <- NULL;
    }
    else{
        if(xregDo=="select"){
            # If this has not happened by chance, then switch to optimisation
            if(!is.null(initialXreg) && (initialType=="optimal")){
                warning("Variables selection does not work with the provided initials for explantory variables. We will drop them.",
                        call.=FALSE);
                initialXreg <- NULL;
                initialXregEstimate <- TRUE;
            }
            else{
                xregDo <- "use";
            }
            if(!is.null(persistenceXreg) && any(persistenceXreg!=0)){
                warning(paste0("We cannot do variables selection with the provided smoothing parameters ",
                               "for explantory variables. We will estimate them instead."),
                        call.=FALSE);
                persistenceXreg <- NULL;
            }
            formulaProvided <- NULL;
        }
    }

    # Use alm() in order to fit the preliminary model for xreg
    if(xregExist){
        xregModel <- vector("list",2);

        # If the initials are not provided, estimate them using ALM.
        if(initialXregEstimate){
            initialXregProvided <- FALSE;
            initialXregEstimate <- TRUE;
            # The function returns an ALM model
            xregInitialiser <- function(Etype,distribution,formulaProvided,otLogical,responseName){
                # Fix the default distribution for ALM
                if(distribution=="default"){
                    distribution <- switch(Etype,
                                           "A"="dnorm",
                                           "M"="dlnorm");
                }
                else if(distribution=="dllaplace"){
                    distribution <- "dlaplace";
                    Etype <- "M";
                }
                else if(distribution=="dls"){
                    distribution <- "ds";
                    Etype <- "M";
                }
                # Return the estimated model based on the provided xreg
                if(is.null(formulaProvided)){
                    if(Etype=="M" && any(distribution==c("dnorm","dlogis","dlaplace","dt","ds","dalaplace"))){
                        formulaProvided <- as.formula(paste0("log(`",responseName,"`)~."));
                    }
                    else{
                        formulaProvided <- as.formula(paste0("`",responseName,"`~."));
                    }
                }
                return(do.call(alm,list(formula=formulaProvided,data=xregData,distribution=distribution,subset=otLogical)))
            }
            # Extract names and form a proper matrix for the regression
            if(!is.null(formulaProvided)){
                formulaProvided <- as.formula(formulaProvided);
                responseName <- all.vars(formulaProvided)[1];
            }

            # If this is not a matrix / data.frame, then convert to one
            if(!is.data.frame(xreg) && !is.matrix(xreg)){
                xreg <- as.data.frame(xreg);
            }

            xregNames <- c(responseName,colnames(xreg));
            if(nrow(xreg)>=obsInSample){
                xregData <- cbind(yInSample,as.data.frame(xreg[1:obsInSample,,drop=FALSE]));
            }
            else{
                stop(paste0("xreg contains less observations than needed: ", nrow(xreg),
                            " instead of ", obsInSample), call.=FALSE);
            }
            colnames(xregData) <- xregNames;

            if(Etype!="Z"){
                testModel <- xregInitialiser(Etype,distribution,formulaProvided,otLogical,responseName);
                if(Etype=="A"){
                    xregModel[[1]]$initialXreg <- testModel$coefficients[-1];
                    xregModel[[1]]$other <- testModel$other;
                }
                else{
                    xregModel[[2]]$initialXreg <- testModel$coefficients[-1];
                    xregModel[[2]]$other <- testModel$other;
                }
            }
            # If we are selecting the appropriate error, produce two models: for "M" and for "A"
            else{
                # Additive model
                testModel <- xregInitialiser("A",distribution,formulaProvided,otLogical,responseName);
                xregModel[[1]]$initialXreg <- testModel$coefficients[-1];
                xregModel[[1]]$other <- testModel$other;
                # Multiplicative model
                testModel[] <- xregInitialiser("M",distribution,formulaProvided,otLogical,responseName);
                xregModel[[2]]$initialXreg <- testModel$coefficients[-1];
                xregModel[[2]]$other <- testModel$other;
            }

            # Write down the number and names of parameters
            xregNumber <- ncol(testModel$data)-1;
            xregData <- testModel$data[,-1,drop=FALSE];
            xregNames <- names(coef(testModel))[-1];
            formulaProvided <- formula(testModel);
            responseName <- formulaProvided[[2]];
            # This is needed in order to succesfully expand the data
            formulaProvided[[2]] <- NULL;

            obsXreg <- nrow(xreg);
            # If there are more xreg values than the obsAll, redo stuff and use them
            if(obsXreg>=obsAll){
                xregData <- model.frame(formulaProvided,data=as.data.frame(xreg));
                xregData <- as.matrix(model.matrix(xregData,data=xregData))[1:obsAll,xregNames,drop=FALSE];
            }
            # If there are less xreg observations than obsAll, use Naive
            else{
                warning(paste0("The xreg has ",obsXreg," observations, while ",obsAll," are needed. ",
                               "Using the last available values as future ones."),
                        call.=FALSE);
                newnRows <- obsAll-obsXreg;
                xregData <- model.frame(formulaProvided,data=as.data.frame(xreg));
                xregData <- as.matrix(model.matrix(xregData,data=xregData))[,xregNames,drop=FALSE];
                xregData <- rbind(xregData,matrix(rep(tail(xregData,1),each=newnRows),newnRows,xregNumber));
            }
            formulaProvided <- formula(testModel);
        }
        else{
            initialXregProvided <- TRUE;

            xregModel[[1]]$initialXreg <- initialXreg;
            if(Etype=="Z"){
                xregModel[[2]]$initialXreg <- initialXreg;
            }

            # Write down the number and names of parameters
            if(nrow(xreg)>obsAll){
                xregData <- xreg[1:obsAll,];
            }
            else if(nrow(xreg)<obsAll){
                stop("The xreg contains less observations than the in-sample. Cannot proceed.",call.=FALSE);
            }
            else{
                xregData <- xreg;
            }
            xregNumber <- ncol(xregData);
            xregNames <- names(xregModel[[1]]$initialXreg);
            parametersNumber[2,2] <- parametersNumber[2,2] + xregNumber;
        }

        # Process the persistence for xreg
        if(!is.null(persistenceXreg)){
            if(length(persistenceXreg)!=xregNumber && length(persistenceXreg)!=1){
                warning("The length of the provided persistence for the xreg variables is wrong. Reverting to the estimation.",
                        call.=FALSE);
                persistenceXreg <- rep(0.5,xregNumber);
                persistenceXregProvided <- FALSE;
                persistenceXregEstimate <- TRUE;
            }
            else if(length(persistenceXreg)==1){
                persistenceXreg <- rep(persistenceXreg,xregNumber);
                persistenceXregProvided <- TRUE;
                persistenceXregEstimate <- FALSE;
            }
            else{
                persistenceXregProvided <- TRUE;
                persistenceXregEstimate <- FALSE;
            }
        }
        else{
            if(xregDo=="adapt"){
                persistenceXreg <- rep(0.05,xregNumber);
                persistenceXregProvided <- FALSE;
                persistenceXregEstimate <- TRUE;
            }
            else{
                persistenceXreg <- rep(0,xregNumber);
                persistenceXregProvided <- FALSE;
                persistenceXregEstimate <- FALSE;
            }
        }
        lagsModelAll <- matrix(c(lagsModelAll,rep(1,xregNumber)),ncol=1);
        # If there's only one explanatory variable, then there's nothing to select
        if(xregNumber==1){
            xregDo[] <- "use";
        }

        # The gsub is needed in order to remove accidental special characters
        colnames(xregData) <- gsub("\`","",colnames(xregData),ignore.case=TRUE);
        xregNames[] <- gsub("\`","",xregNames,ignore.case=TRUE);
    }
    else{
        initialXregProvided <- FALSE;
        initialXregEstimate <- FALSE;
        persistenceXregProvided <- FALSE;
        persistenceXregEstimate <- FALSE;
        xregModel <- NULL;
        xregData <- NULL;
        xregNumber <- 0;
        xregNames <- NULL;
        if(is.null(formulaProvided)){
            if(Etype=="M" && any(distribution==c("dnorm","dlogis","dlaplace","dt","ds","dalaplace"))){
                formulaProvided <- as.formula(paste0("log(`",responseName,"`)~."));
            }
            else{
                formulaProvided <- as.formula(paste0("`",responseName,"`~."));
            }
        }
    }
    # Remove xreg, just to preserve some memory
    rm(xreg);

    #### Conclusions about the initials ####
    # Make sure that only important elements are estimated.
    if(Ttype=="N"){
        initialTrendEstimate <- FALSE;
        initialTrend <- NULL;
    }
    if(Stype=="N"){
        initialSeasonalEstimate <- FALSE;
        initialSeasonal <- NULL;
    }
    if(!arimaModel){
        initialArimaEstimate <- FALSE;
        initialArima <- NULL;
    }
    if(!xregExist){
        initialXregEstimate <- FALSE;
        initialXreg <- NULL;
    }

    # If we don't need to estimate anything, flag initialEstimate
    if(!any(c(initialLevelEstimate, (initialTrendEstimate & modelIsTrendy),
              (initialSeasonalEstimate & Stype!="N"),
              (initialArimaEstimate & arimaModel),
              (initialXregEstimate & xregExist)))){
        initialEstimate[] <- FALSE;
    }
    else{
        initialEstimate[] <- TRUE;
    }
    # If at least something is provided, flag it as "provided"
    if(!all(c(initialLevelEstimate, (initialTrendEstimate & modelIsTrendy),
              (initialSeasonalEstimate & Stype!="N"),
              (initialArimaEstimate & arimaModel),
              (initialXregEstimate & xregExist)))){
        initialType[] <- "provided";
    }

    # Observations in the states matrix
    # Define the number of cols that should be in the matvt
    obsStates <- obsInSample + lagsModelMax*switch(initialType,
                                                   "backcasting"=2,
                                                   1);


    # Check if multiplicative models can be fitted
    allowMultiplicative <- !((any(yInSample<=0) && !occurrenceModel) || (occurrenceModel && any(yInSample<0)));

    # Clean the pool of models if only additive are allowed
    if(!allowMultiplicative && !is.null(modelsPool)){
        modelsPoolMultiplicative <- ((substr(modelsPool,1,1)=="M") |
                                         substr(modelsPool,2,2)=="M" |
                                         substr(modelsPool,nchar(modelsPool),nchar(modelsPool))=="M");
        if(any(modelsPoolMultiplicative)){
            modelsPool <- modelsPool[!modelsPoolMultiplicative];

            if(!any(model==c("PPP","FFF"))){
                warning("Only additive models are allowed for your data. Amending the pool.",
                        call.=FALSE);
            }
        }
    }
    if(any(model==c("PPP","FFF"))){
        model <- "ZZZ";
    }

    # Update the number of parameters
    if(occurrenceModelProvided){
        parametersNumber[2,3] <- nparam(oesModel);
        pForecast <- c(forecast(oesModel, h=h)$mean);
    }

    #### Information Criteria ####
    ic <- match.arg(ic,c("AICc","AIC","BIC","BICc"));
    ICFunction <- switch(ic,
                         "AIC"=AIC,
                         "AICc"=AICc,
                         "BIC"=BIC,
                         "BICc"=BICc);

    #### Bounds for the smoothing parameters ####
    bounds <- match.arg(bounds,c("usual","admissible","none"));


    #### Checks for the potential number of degrees of freedom ####
    # This is needed in order to make the function work on small samples
    # scale parameter, smoothing parameters and phi
    nParamMax <- (1 + persistenceLevelEstimate + persistenceTrendEstimate*modelIsTrendy +
                      sum(persistenceSeasonalEstimate)*modelIsSeasonal +
                      phiEstimate +
                      # Number of ETS initials
                      (sum(lagsModelAll)-xregNumber-initialArimaNumber)*(initialType=="optimal") +
                      # ARIMA components: initials + parameters
                      arimaModel*(initialArimaNumber*(initialType=="optimal") + sum(arOrders) + sum(maOrders)) +
                      # Xreg initials and smoothing parameters
                      xregNumber*(initialXregEstimate+persistenceXregEstimate));

    # If the sample is smaller than the number of parameters
    if(obsNonzero <= nParamMax){
        # If there is ARIMA terms, remove them
        if(arimaModel){
            warning("We don't have enough observations to fit ETS with ARIMA terms. We will construct the simple ETS.",
                    call.=FALSE);
            arRequired <- iRequired <- maRequired <- arimaModel <- FALSE;
            arOrders <- iOrders <- maOrders <- NULL;
            nonZeroARI <- nonZeroMA <- lagsModelARIMA <- NULL;
            componentsNamesARIMA <- NULL;
            initialArimaNumber <- componentsNumberARIMA <- 0;
            lagsModelAll <- lagsModelAll[-c(componentsNumberETS+c(1:componentsNumberARIMA)),,drop=FALSE];
            lagsModelMax <- max(lagsModelAll);

            nParamMax[] <- (1 + persistenceLevelEstimate + persistenceTrendEstimate*modelIsTrendy +
                                sum(persistenceSeasonalEstimate)*modelIsSeasonal +
                                phiEstimate +
                                # Number of ETS initials
                                (sum(lagsModelAll)-xregNumber)*(initialType=="optimal") +
                                # Xreg initials and smoothing parameters
                                xregNumber*(initialXregEstimate+persistenceXregEstimate));
        }
    }

    # If the sample is still smaller than the number of parameters
    if(obsNonzero <= nParamMax){
        nParamExo <- xregNumber*(initialXregEstimate+persistenceXregEstimate);
        if(!silent){
            message(paste0("Number of non-zero observations is ",obsNonzero,
                           ", while the maximum number of parameters to estimate is ", nParamMax,".\n",
                           "Updating pool of models."));
        }

        # If the number of observations is still enough for the model selection and the pool is not specified
        if(obsNonzero > (3 + nParamExo) && is.null(modelsPool) && any(modelDo==c("select","combine"))){
            # We have enough observations for local level model
            modelsPool <- c("ANN");
            if(allowMultiplicative){
                modelsPool <- c(modelsPool,"MNN");
            }
            # We have enough observations for trend model
            if(obsNonzero > (5 + nParamExo)){
                modelsPool <- c(modelsPool,"AAN");
                if(allowMultiplicative){
                    modelsPool <- c(modelsPool,"AMN","MAN","MMN");
                }
            }
            # We have enough observations for damped trend model
            if(obsNonzero > (6 + nParamExo)){
                modelsPool <- c(modelsPool,"AAdN");
                if(allowMultiplicative){
                    modelsPool <- c(modelsPool,"AMdN","MAdN","MMdN");
                }
            }
            # We have enough observations for seasonal model
            if((obsNonzero > (2*lagsModelMax)) && lagsModelMax!=1){
                modelsPool <- c(modelsPool,"ANA");
                if(allowMultiplicative){
                    modelsPool <- c(modelsPool,"ANM","MNA","MNM");
                }
            }
            # We have enough observations for seasonal model with trend
            if((obsNonzero > (6 + lagsModelMax + nParamExo)) &&
               (obsNonzero > 2*lagsModelMax) && lagsModelMax!=1){
                modelsPool <- c(modelsPool,"AAA");
                if(allowMultiplicative){
                    modelsPool <- c(modelsPool,"AAM","AMA","AMM","MAA","MAM","MMA","MMM");
                }
            }

            warning("Not enough of non-zero observations for the fit of ETS(",model,")! Fitting what we can...",
                    call.=FALSE);
            if(modelDo=="combine"){
                model <- "CNN";
                if(length(modelsPool)>2){
                    model <- "CCN";
                }
                if(length(modelsPool)>10){
                    model <- "CCC";
                }
            }
            else{
                modelDo <- "select"
                model <- "ZZZ";
            }
        }
        # If the pool is provided (so, select / combine), amend it
        else if(obsNonzero > (3 + nParamExo) && !is.null(modelsPool)){
            # We don't have enough observations for seasonal models with damped trend
            if((obsNonzero <= (6 + lagsModelMax + 1 + nParamExo))){
                modelsPool <- modelsPool[!(nchar(modelsPool)==4 &
                                               substr(modelsPool,nchar(modelsPool),nchar(modelsPool))=="A")];
                modelsPool <- modelsPool[!(nchar(modelsPool)==4 &
                                               substr(modelsPool,nchar(modelsPool),nchar(modelsPool))=="M")];
            }
            # We don't have enough observations for seasonal models with trend
            if((obsNonzero <= (5 + lagsModelMax + 1 + nParamExo))){
                modelsPool <- modelsPool[!(substr(modelsPool,2,2)!="N" &
                                               substr(modelsPool,nchar(modelsPool),nchar(modelsPool))!="N")];
            }
            # We don't have enough observations for seasonal models
            if(obsNonzero <= 2*lagsModelMax){
                modelsPool <- modelsPool[substr(modelsPool,nchar(modelsPool),nchar(modelsPool))=="N"];
            }
            # We don't have enough observations for damped trend
            if(obsNonzero <= (6 + nParamExo)){
                modelsPool <- modelsPool[nchar(modelsPool)!=4];
            }
            # We don't have enough observations for any trend
            if(obsNonzero <= (5 + nParamExo)){
                modelsPool <- modelsPool[substr(modelsPool,2,2)=="N"];
            }

            modelsPool <- unique(modelsPool);
            warning("Not enough of non-zero observations for the fit of ETS(",model,")! Fitting what we can...",
                    call.=FALSE);
            if(modelDo=="combine"){
                model <- "CNN";
                if(length(modelsPool)>2){
                    model <- "CCN";
                }
                if(length(modelsPool)>10){
                    model <- "CCC";
                }
            }
            else{
                modelDo <- "select"
                model <- "ZZZ";
            }
        }
        # If the model needs to be estimated / used, not selected
        else if(obsNonzero > (3 + nParamExo) && any(modelDo==c("estimate","use"))){
            # We don't have enough observations for seasonal models with damped trend
            if((obsNonzero <= (6 + lagsModelMax + 1 + nParamExo))){
                model <- model[!(nchar(model)==4 &
                                               substr(model,nchar(model),nchar(model))=="A")];
                model <- model[!(nchar(model)==4 &
                                               substr(model,nchar(model),nchar(model))=="M")];
            }
            # We don't have enough observations for seasonal models with trend
            if((obsNonzero <= (5 + lagsModelMax + 1 + nParamExo))){
                model <- model[!(substr(model,2,2)!="N" &
                                               substr(model,nchar(model),nchar(model))!="N")];
            }
            # We don't have enough observations for seasonal models
            if(obsNonzero <= 2*lagsModelMax){
                model <- model[substr(model,nchar(model),nchar(model))=="N"];
            }
            # We don't have enough observations for damped trend
            if(obsNonzero <= (6 + nParamExo)){
                model <- model[nchar(model)!=4];
            }
            # We don't have enough observations for any trend
            if(obsNonzero <= (5 + nParamExo)){
                model <- model[substr(model,2,2)=="N"];
            }
        }
        # Extreme cases of small samples
        else if(obsNonzero==4){
            if(any(Etype==c("A","M"))){
                modelDo <- "estimate";
                Ttype <- "N";
                Stype <- "N";
            }
            else{
                modelsPool <- c("ANN");
                if(allowMultiplicative){
                    modelsPool <- c(modelsPool,"MNN");
                }
                modelDo <- "select";
                model <- "ZZZ";
                Etype <- "Z";
                Ttype <- "N";
                Stype <- "N";
                warning("You have a very small sample. The only available model is level model.",
                        call.=FALSE);
            }
            phiEstimate <- FALSE;
        }
        # Even smaller sample
        else if(obsNonzero==3){
            if(any(Etype==c("A","M"))){
                modelDo <- "estimate";
                Ttype <- "N";
                Stype <- "N";
                model <- paste0(Etype,"NN");
            }
            else{
                modelsPool <- c("ANN");
                if(allowMultiplicative){
                    modelsPool <- c(modelsPool,"MNN");
                }
                modelDo <- "select";
                model <- "ZNN";
                Etype <- "Z";
                Ttype <- "N";
                Stype <- "N";
            }
            persistence <- 0;
            names(persistence) <- "level";
            persistenceEstimate <- persistenceLevelEstimate <- FALSE;
            warning("We did not have enough of non-zero observations, so persistence value was set to zero.",
                    call.=FALSE);
            phiEstimate <- FALSE;
        }
        # Can it be even smaller?
        else if(obsNonzero==2){
            modelsPool <- NULL;
            persistence <- 0;
            names(persistence) <- "level";
            persistenceEstimate <- persistenceLevelEstimate <- FALSE;
            initialLevel <- mean(yInSample);
            initialType <- "provided";
            initialEstimate <- initialLevelEstimate <- FALSE;
            warning("We did not have enough of non-zero observations, so persistence value was set to zero and initial was preset.",
                    call.=FALSE);
            modelDo <- "use";
            model <- "ANN";
            Etype <- "A";
            Ttype <- "N";
            Stype <- "N";
            phiEstimate <- FALSE;
            parametersNumber[1,1] <- 0;
            parametersNumber[2,1] <- 2;
        }
        # And how about now?!
        else if(obsNonzero==1){
            modelsPool <- NULL;
            persistence <- 0;
            names(persistence) <- "level";
            persistenceEstimate <- persistenceLevelEstimate <- FALSE;
            initialLevel <- yInSample[yInSample!=0];
            initialType <- "provided";
            initialEstimate <- initialLevelEstimate <- FALSE;
            warning("We did not have enough of non-zero observations, so we used Naive.",call.=FALSE);
            modelDo <- "nothing"
            model <- "ANN";
            Etype <- "A";
            Ttype <- "N";
            Stype <- "N";
            phiEstimate <- FALSE;
            parametersNumber[1,1] <- 0;
            parametersNumber[2,1] <- 2;
        }
        # Only zeroes in the data...
        else if(obsNonzero==0 && obsInSample>1){
            modelsPool <- NULL;
            persistence <- 0;
            names(persistence) <- "level";
            persistenceEstimate <- persistenceLevelEstimate <- FALSE;
            initialLevel <- 0;
            initialType <- "provided";
            initialEstimate <- initialLevelEstimate <- FALSE;
            occurrenceModelProvided <- occurrenceModel <- FALSE;
            occurrence <- "none";
            warning("You have a sample with zeroes only. Your forecast will be zero.",call.=FALSE);
            modelDo <- "nothing"
            model <- "ANN";
            Etype <- "A";
            Ttype <- "N";
            Stype <- "N";
            phiEstimate <- FALSE;
            parametersNumber[1,1] <- 0;
            parametersNumber[2,1] <- 2;
        }
        # If you don't have observations, then fuck off!
        else{
            stop("Not enough observations... Even for fitting of ETS('ANN')!",call.=FALSE);
        }
    }
    # Reset the maximum lag. This is in order to take potential changes into account
    lagsModelMax[] <- max(lagsModelAll);

    #### Process ellipsis ####
    # Parameters for the optimiser
    if(is.null(ellipsis$maxeval)){
        if(arimaModel){
            maxeval <- 1000;
        }
        else{
            maxeval <- 200;
        }
        # This is heuristic. If you have higher seasonal lags, use more iterations.
        if(lagsModelMax>12){
            maxeval[] <- maxeval/20 * lagsModelMax;
        }
    }
    else{
        maxeval <- ellipsis$maxeval;
    }
    if(is.null(ellipsis$maxtime)){
        maxtime <- -1;
    }
    else{
        maxtime <- ellipsis$maxtime;
    }
    if(is.null(ellipsis$xtol_rel)){
        xtol_rel <- 1E-6;
    }
    else{
        xtol_rel <- ellipsis$xtol_rel;
    }
    if(is.null(ellipsis$xtol_abs)){
        xtol_abs <- 0;
    }
    else{
        xtol_abs <- ellipsis$xtol_abs;
    }
    if(is.null(ellipsis$algorithm)){
        algorithm <- "NLOPT_LN_NELDERMEAD";
    }
    else{
        algorithm <- ellipsis$algorithm;
    }
    if(is.null(ellipsis$print_level)){
        print_level <- 0;
    }
    else{
        print_level <- ellipsis$print_level;
    }
    # The following three arguments are used for the function itself, not the options
    if(is.null(ellipsis$lb)){
        lb <- NULL;
    }
    else{
        lb <- ellipsis$lb;
    }
    if(is.null(ellipsis$ub)){
        ub <- NULL;
    }
    else{
        ub <- ellipsis$ub;
    }
    if(is.null(ellipsis$B)){
        B <- NULL;
    }
    else{
        B <- ellipsis$B;
    }
    # Additional parameter for dalaplace, LASSO and dt
    if(is.null(ellipsis$lambda)){
        if(loss=="likelihood" && any(distribution==c("dt","dalaplace"))){
            lambdaEstimate <- TRUE;
        }
        else{
            lambdaEstimate <- FALSE;
        }
        lambda <- 0.5;
    }
    else{
        lambdaEstimate <- FALSE;
        lambda <- ellipsis$lambda;
    }
    # Fisher Information
    if(is.null(ellipsis$FI)){
        FI <- FALSE;
    }
    else{
        FI <- ellipsis$FI;
    }

    # See if the estimation of the model is not needed
    if(!any(persistenceEstimate,phiEstimate, (initialType!="backcasting")&initialEstimate,
            arimaModel, lambdaEstimate)){
        modelDo <- "use";
    }

    #### Return the values to the previous environment ####
    ### Actuals
    assign("y",y,ParentEnvironment);
    assign("yHoldout",yHoldout,ParentEnvironment);
    assign("yInSample",yInSample,ParentEnvironment);
    assign("yNAValues",yNAValues,ParentEnvironment);

    ### Index and all related structure variables
    assign("yClasses",yClasses,ParentEnvironment);
    assign("yIndex",yIndex,ParentEnvironment);
    assign("yInSampleIndex",yInSampleIndex,ParentEnvironment);
    assign("yForecastIndex",yForecastIndex,ParentEnvironment);
    assign("yFrequency",yFrequency,ParentEnvironment);
    assign("yStart",yStart,ParentEnvironment);
    assign("yForecastStart",yForecastStart,ParentEnvironment);

    # The rename of the variable is needed for the hessian to work
    assign("horizon",h,ParentEnvironment);
    assign("h",h,ParentEnvironment);
    assign("holdout",holdout,ParentEnvironment);

    ### Number of observations and parameters
    assign("obsInSample",obsInSample,ParentEnvironment);
    assign("obsAll",obsAll,ParentEnvironment);
    assign("obsStates",obsStates,ParentEnvironment);
    assign("obsNonzero",obsNonzero,ParentEnvironment);
    assign("obsZero",obsZero,ParentEnvironment);
    assign("parametersNumber",parametersNumber,ParentEnvironment);

    ### Model type
    assign("model",model,ParentEnvironment);
    assign("Etype",Etype,ParentEnvironment);
    assign("Ttype",Ttype,ParentEnvironment);
    assign("Stype",Stype,ParentEnvironment);
    assign("modelIsTrendy",modelIsTrendy,ParentEnvironment);
    assign("modelIsSeasonal",modelIsSeasonal,ParentEnvironment);
    assign("modelsPool",modelsPool,ParentEnvironment);
    assign("damped",damped,ParentEnvironment);
    assign("modelDo",modelDo,ParentEnvironment);
    assign("allowMultiplicative",allowMultiplicative,ParentEnvironment);

    ### Numbers and names of components
    assign("componentsNumberETS",componentsNumberETS,ParentEnvironment);
    assign("componentsNamesETS",componentsNamesETS,ParentEnvironment);
    assign("componentsNumberETSNonSeasonal",componentsNumberETS-componentsNumberETSSeasonal,ParentEnvironment);
    assign("componentsNumberETSSeasonal",componentsNumberETSSeasonal,ParentEnvironment);
    # The number and names of ARIMA components
    assign("componentsNumberARIMA",componentsNumberARIMA,ParentEnvironment);
    assign("componentsNamesARIMA",componentsNamesARIMA,ParentEnvironment);

    ### Lags
    # This is the original vector of lags, modified for the level components.
    # This can be used in ARIMA
    assign("lags",lags,ParentEnvironment);
    # This is the vector of lags of ETS components
    assign("lagsModel",lagsModel,ParentEnvironment);
    # This is the vector of seasonal lags
    assign("lagsModelSeasonal",lagsModelSeasonal,ParentEnvironment);
    # This is the vector of lags for ARIMA components (not lags of ARIMA)
    assign("lagsModelARIMA",lagsModelARIMA,ParentEnvironment);
    # This is the vector of all the lags of model (ETS + ARIMA + X)
    assign("lagsModelAll",lagsModelAll,ParentEnvironment);
    # This is the maximum lag
    assign("lagsModelMax",lagsModelMax,ParentEnvironment);

    ### Persistence
    assign("persistence",persistence,ParentEnvironment);
    assign("persistenceEstimate",persistenceEstimate,ParentEnvironment);
    assign("persistenceLevel",persistenceLevel,ParentEnvironment);
    assign("persistenceLevelEstimate",persistenceLevelEstimate,ParentEnvironment);
    assign("persistenceTrend",persistenceTrend,ParentEnvironment);
    assign("persistenceTrendEstimate",persistenceTrendEstimate,ParentEnvironment);
    assign("persistenceSeasonal",persistenceSeasonal,ParentEnvironment);
    assign("persistenceSeasonalEstimate",persistenceSeasonalEstimate,ParentEnvironment);
    assign("persistenceXreg",persistenceXreg,ParentEnvironment);
    assign("persistenceXregEstimate",persistenceXregEstimate,ParentEnvironment);
    assign("persistenceXregProvided",persistenceXregProvided,ParentEnvironment);

    ### phi
    assign("phi",phi,ParentEnvironment);
    assign("phiEstimate",phiEstimate,ParentEnvironment);

    ### Initials
    assign("initial",initial,ParentEnvironment);
    assign("initialType",initialType,ParentEnvironment);
    assign("initialEstimate",initialEstimate,ParentEnvironment);
    assign("initialLevel",initialLevel,ParentEnvironment);
    assign("initialLevelEstimate",initialLevelEstimate,ParentEnvironment);
    assign("initialTrend",initialTrend,ParentEnvironment);
    assign("initialTrendEstimate",initialTrendEstimate,ParentEnvironment);
    assign("initialSeasonal",initialSeasonal,ParentEnvironment);
    assign("initialSeasonalEstimate",initialSeasonalEstimate,ParentEnvironment);
    assign("initialArima",initialArima,ParentEnvironment);
    assign("initialArimaEstimate",initialArimaEstimate,ParentEnvironment);
    # Number of initials that the ARIMA has (either provided or to estimate)
    assign("initialArimaNumber",initialArimaNumber,ParentEnvironment);
    assign("initialXreg",initialXreg,ParentEnvironment);
    assign("initialXregEstimate",initialXregEstimate,ParentEnvironment);
    assign("initialXregProvided",initialXregProvided,ParentEnvironment);

    ### Occurrence model
    assign("oesModel",oesModel,ParentEnvironment);
    assign("occurrenceModel",occurrenceModel,ParentEnvironment);
    assign("occurrenceModelProvided",occurrenceModelProvided,ParentEnvironment);
    assign("occurrence",occurrence,ParentEnvironment);
    assign("pFitted",pFitted,ParentEnvironment);
    assign("pForecast",pForecast,ParentEnvironment);
    assign("ot",ot,ParentEnvironment);
    assign("otLogical",otLogical,ParentEnvironment);

    ### Distribution, loss, bounds and IC
    assign("distribution",distribution,ParentEnvironment);
    assign("loss",loss,ParentEnvironment);
    assign("lossFunction",lossFunction,ParentEnvironment);
    assign("multisteps",multisteps,ParentEnvironment);
    assign("ic",ic,ParentEnvironment);
    assign("ICFunction",ICFunction,ParentEnvironment);
    assign("bounds",bounds,ParentEnvironment);

    ### ARIMA components
    assign("arOrders",arOrders,ParentEnvironment);
    assign("iOrders",iOrders,ParentEnvironment);
    assign("maOrders",maOrders,ParentEnvironment);
    assign("arimaModel",arimaModel,ParentEnvironment);
    assign("arRequired",arRequired,ParentEnvironment);
    assign("iRequired",iRequired,ParentEnvironment);
    assign("maRequired",maRequired,ParentEnvironment);
    assign("arEstimate",arEstimate,ParentEnvironment);
    assign("maEstimate",maEstimate,ParentEnvironment);
    assign("nonZeroARI",nonZeroARI,ParentEnvironment);
    assign("nonZeroMA",nonZeroMA,ParentEnvironment);

    ### Explanatory variables
    assign("xregDo",xregDo,ParentEnvironment);
    assign("xregExist",xregExist,ParentEnvironment);
    assign("xregModel",xregModel,ParentEnvironment);
    assign("xregData",xregData,ParentEnvironment);
    assign("xregNumber",xregNumber,ParentEnvironment);
    assign("xregNames",xregNames,ParentEnvironment);
    assign("formula",formulaProvided,ParentEnvironment);

    ### Ellipsis thingies
    # Optimisation related
    assign("maxeval",maxeval,ParentEnvironment);
    assign("maxtime",maxtime,ParentEnvironment);
    assign("xtol_rel",xtol_rel,ParentEnvironment);
    assign("xtol_abs",xtol_abs,ParentEnvironment);
    assign("algorithm",algorithm,ParentEnvironment);
    assign("print_level",print_level,ParentEnvironment);
    assign("B",B,ParentEnvironment);
    assign("lb",lb,ParentEnvironment);
    assign("ub",ub,ParentEnvironment);
    # Additional parameters
    assign("lambda",lambda,ParentEnvironment);
    assign("lambdaEstimate",lambdaEstimate,ParentEnvironment);
    assign("FI",FI,ParentEnvironment);
}
