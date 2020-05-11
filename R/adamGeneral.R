parametersChecker <- function(y, model, lags, formulaProvided, orders,
                              persistence, phi, initial,
                              distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                                             "dlnorm","dllaplace","dls","dinvgauss"),
                              loss, h, holdout,occurrence,
                              ic=c("AICc","AIC","BIC","BICc"), bounds=c("traditional","admissible","none"),
                              xreg, xregDo, xregInitial, xregPersistence, responseName,
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
            yIndex <- yIndex[1] + c(1:length(y[[1]])) * diff(yIndex)[1];
        }
    }
    else{
        yIndex <- time(y);
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
        else{
            if(ncol(y)>1){
                xreg <- y[,-1];
            }
            y <- y[,1];
        }
    }

    # Substitute NAs with mean values.
    if(any(is.na(y))){
        warning("Data contains NAs. The values will be ignored during the model construction.",call.=FALSE);
        yNAValues <- is.na(y);
        y[yNAValues] <- na.interp(y)[yNAValues];
    }
    else{
        yNAValues <- is.na(y);
    }

    # Define obs, the number of observations of in-sample
    obsAll <- length(y) + (1 - holdout)*h;
    obsInSample <- length(y) - holdout*h;

    # If this is just a numeric variable, use ts class
    if(all(yClasses=="integer") || all(yClasses=="data.frame") || all(yClasses=="matrix")){
        yClasses <- "ts";
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
        yForecastStart <- yIndex[obsInSample]+diff(yIndex)[1];
        yInSampleIndex <- yIndex;
        yForecastIndex <- yIndex[obsInSample]+diff(yIndex)[1]*c(1:max(h,1));
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
        Etype <- substring(model,1,1);
        Ttype <- substring(model,2,2);
        Stype <- substring(model,4,4);
        damped <- TRUE;
        if(substring(model,3,3)!="d"){
            message(paste0("You have defined a strange model: ",model));
            model <- paste0(Etype,Ttype,"d",Stype);
        }
    }
    else if(nchar(model)==3){
        Etype <- substring(model,1,1);
        Ttype <- substring(model,2,2);
        Stype <- substring(model,3,3);
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

    #### Check the components of model ####
    componentsNames <- "level";
    componentsNumber <- 1;
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
    if(Ttype!="N"){
        componentsNames <- c(componentsNames,"trend");
        componentsNumber[] <- componentsNumber+1;
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
        if(sum(arOrders)>0){
            arRequired <- TRUE;
        }
        else{
            arRequired <- FALSE;
        }

        # See if I is needed
        if(sum(iOrders)>0){
            iRequired <- TRUE;
        }
        else{
            iRequired <- FALSE;
        }

        # See if I is needed
        if(sum(maOrders)>0){
            maRequired <- TRUE;
        }
        else{
            maRequired <- FALSE;
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

        # If zeroes are defined for some orders, drop them.
        # if(any((arOrders + iOrders + maOrders)==0)){
        #     orders2leave <- (arOrders + iOrders + maOrders)!=0;
        #     if(all(!orders2leave)){
        #         orders2leave <- lags==min(lags);
        #     }
        #     arOrders <- arOrders[orders2leave];
        #     iOrders <- iOrders[orders2leave];
        #     maOrders <- maOrders[orders2leave];
        #     lags <- lags[orders2leave];
        # }
    }
    else{
        arimaModel <- FALSE;
        arRequired <- arEstimate <- FALSE;
        iRequired <- FALSE;
        maRequired <- maEstimate <- FALSE;
        lagsModelARIMA <- initialNumberARIMA <- 0;
        componentsNumberARIMA <- 0;
        componentsNamesARIMA <- NULL;
    }

    if(arimaModel){
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
        lagsModelARIMA <- matrix(sort(unique(c(nonZeroARI,nonZeroMA))),ncol=1);
        nonZeroARI <- cbind(nonZeroARI,which(lagsModelARIMA %in% nonZeroARI)-1);
        nonZeroMA <- cbind(nonZeroMA,which(lagsModelARIMA %in% nonZeroMA)-1);

        componentsNumberARIMA <- length(lagsModelARIMA);
        componentsNamesARIMA <- paste0("State",c(1:componentsNumberARIMA));
        initialNumberARIMA <- sum(lagsModelARIMA)

        if(obsInSample < componentsNumberARIMA){
            warning(paste0("In-sample size is ",obsInSample,", while number of ARIMA components is ",componentsNumberARIMA,
                           ". Cannot fit the model."),call.=FALSE)
            stop("Not enough observations for such a complicated model.",call.=FALSE);
        }

        # Check the provided parameters for AR and MA

        # Check the provided initials

    }
    else{
        componentsNumberARIMA <- 0;
        nonZeroARI <- NULL;
        nonZeroMA <- NULL;
        arOrders <- NULL;
        iOrders <- NULL;
        maOrders <- NULL;
    }

    # If we have a trend add one more lag
    if(Ttype!="N"){
        lags <- c(1,lags);
    }
    # If we don't have seasonality, remove seasonal lag
    if(Stype=="N" & any(lags>1)){
        lags <- lags[lags==1];
    }

    # Lags of the model
    lagsModel <- matrix(lags,ncol=1);
    lagsModelMax <- max(lagsModel);
    lagsLength <- length(lags);

    modelIsSeasonal <- Stype!="N";
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
                           "Switching to non-seasonal model: ETS(",substring(model,1,nchar(model)-1),"N)"));
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
        componentsNames <- c(componentsNames,"seasonal");
        componentsNumber[] <- componentsNumber+1;
        componentsNumberSeasonal <- 1;
    }
    else{
        componentsNumberSeasonal <- 0;
    }

    # Check, whether the number of lags and the number of components are the same
    if(lagsLength>componentsNumber){
        if(Stype!="N"){
            componentsNames <- c(componentsNames[-length(componentsNames)],paste0("seasonal",c(1:(lagsLength-componentsNumber-1))));
            componentsNumberSeasonal[] <- lagsLength-componentsNumber+1;
            componentsNumber[] <- lagsLength;
        }
        else{
            lagsModel <- matrix(lags[1:componentsNumber],ncol=1);
            lagsModelMax <- max(lagsModel);
            lagsLength <- length(lagsModel);
        }
    }
    else if(lagsLength<componentsNumber){
        stop("The number of components of the model is smaller than the number of provided lags", call.=FALSE);
    }

    if(!fast){
        #### Distribution selected ####
        distribution <- match.arg(distribution);
    }

    #### Loss function type ####
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

    #### Persistence provided ####
    if(!is.null(persistence)){
        if((!is.numeric(persistence) | !is.vector(persistence)) & !is.matrix(persistence)){
            warning(paste0("Persistence is not a numeric vector!\n",
                           "Changing to estimation of persistence vector values."),call.=FALSE);
            persistence <- NULL;
            persistenceEstimate <- TRUE;
        }
        else{
            if(modelDo!="estimate"){
                warning(paste0("Predefined persistence vector can only be used with ",
                               "preselected ETS model.\n",
                               "Changing to estimation of persistence vector values."),call.=FALSE);
                persistence <- NULL;
                persistenceEstimate <- TRUE;
            }
            else{
                if(length(persistence)!=lagsLength){
                    warning(paste0("Length of persistence vector is wrong! ",
                                   "It should be ",lagsLength,".\n",
                                   "Changing to estimation of persistence vector values."),
                            call.=FALSE);
                    persistence <- NULL;
                    persistenceEstimate <- TRUE;
                }
                else{
                    persistence <- as.vector(persistence);
                    names(persistence) <- componentsNames;
                    persistenceEstimate <- FALSE;
                    parametersNumber[2,1] <- parametersNumber[2,1] + length(persistence);
                    bounds <- "n";
                }
            }
        }
    }
    else{
        persistenceEstimate <- TRUE;
    }

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

    #### Vector of initial values ####
    # initial type can be: "o" - optimal, "b" - backcasting, "p" - provided.
    if(any(is.character(initial))){
        initialType <- match.arg(initial, c("optimal","backcasting"));
        initialValue <- NULL;
    }
    else if(is.null(initial)){
        if(!silent){
            message("Initial value is not selected. Switching to optimal.");
        }
        initialType <- "optimal";
        initialValue <- NULL;
    }
    else if(!is.null(initial)){
        if(!is.numeric(initial)){
            warning(paste0("Initial vector is not numeric!\n",
                           "Values of initial vector will be estimated."),call.=FALSE);
            initialValue <- NULL;
            initialType <- "optimal";
        }
        else{
            if(modelDo!="estimate"){
                warning(paste0("Predefined initials vector can only be used with preselected ETS model.\n",
                               "Changing to estimation of initials."),call.=FALSE);
                initialValue <- NULL;
                initialType <- "optimal";
            }
            else{
                if(length(initial)!=sum(lags)){
                    warning(paste0("Wrong length of the initial vector. Should be ",sum(lags),
                                   " instead of ",length(initial),".\n",
                                   "Values of initial vector will be estimated."),call.=FALSE);
                    initialValue <- NULL;
                    initialType <- "optimal";
                }
                else{
                    initialType <- "provided";
                    initialValue <- initial;
                    parametersNumber[2,1] <- parametersNumber[2,1] + sum(lags);
                }
            }
        }
    }
    initialEstimate <- initialType=="optimal";

    # Observations in the states matrix
    # Define the number of cols that should be in the matvt
    obsStates <- obsInSample + lagsModelMax*switch(initialType,
                                                   "backcasting"=2,
                                                   1);

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

    #### Lags for ARIMA ####
    if(arimaModel){
        lagsModelAll <- rbind(lagsModel,lagsModelARIMA);
        lagsModelMax <- max(lagsModel);
    }
    else{
        lagsModelAll <- lagsModel;
    }

    #### Explanatory variables: xreg, xregDo, xregInitial, xregPersistence ####
    xregDo <- match.arg(xregDo,c("use","select"));
    xregExist <- !is.null(xreg);
    if(!xregExist){
        xregDo[] <- "use";
        formulaProvided <- NULL;
    }
    else{
        if(xregDo=="select"){
            if(!is.null(xregInitial)){
                warning("Variables selection does not work with the provided initials for explantory variables. We will drop them.",
                        call.=FALSE);
                xregInitial <- NULL;
            }
            if(!is.null(xregPersistence) && any(xregPersistence!=0)){
                warning(paste0("We cannot do variables selection with the provided smoothing parameters ",
                               "for explantory variables. We will estimate them instead."),
                        call.=FALSE);
                xregPersistence <- NULL;
            }
            formulaProvided <- NULL;
        }
    }

    # Use alm() in order to fit the preliminary model for xreg
    if(xregExist){
        xregModel <- vector("list",2);

        # If the initials are not provided, estimate them using ALM.
        if(is.null(xregInitial)){
            xregInitialsProvided <- FALSE;
            xregInitialsEstimate <- TRUE;
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
            xregData <- cbind(yInSample,xreg[1:obsInSample,,drop=FALSE]);
            colnames(xregData) <- xregNames;

            if(Etype!="Z"){
                testModel <- xregInitialiser(Etype,distribution,formulaProvided,otLogical,responseName);
                if(Etype=="A"){
                    xregModel[[1]]$xregInitial <- testModel$coefficients[-1];
                    xregModel[[1]]$other <- testModel$other;
                }
                else{
                    xregModel[[2]]$xregInitial <- testModel$coefficients[-1];
                    xregModel[[2]]$other <- testModel$other;
                }
            }
            # If we are selecting the appropriate error, produce two models: for "M" and for "A"
            else{
                # Additive model
                testModel <- xregInitialiser("A",distribution,formulaProvided,otLogical,responseName);
                xregModel[[1]]$xregInitial <- testModel$coefficients[-1];
                xregModel[[1]]$other <- testModel$other;
                # Multiplicative model
                testModel[] <- xregInitialiser("M",distribution,formulaProvided,otLogical,responseName);
                xregModel[[2]]$xregInitial <- testModel$coefficients[-1];
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
            xregInitialsProvided <- TRUE;
            xregInitialsEstimate <- FALSE;

            xregModel[[1]]$xregInitial <- xregInitial;
            if(Etype=="Z"){
                xregModel[[2]]$xregInitial <- xregInitial;
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
            xregNames <- names(xregModel[[1]]$xregInitial);
        }

        # Process the persistence for xreg
        if(!is.null(xregPersistence)){
            if(length(xregPersistence)!=xregNumber && length(xregPersistence)!=1){
                warning("The length of the provided xregPersistence variables is wrong. Reverting to the estimation.",
                        call.=FALSE);
                xregPersistence <- rep(0.5,xregNumber);
                xregPersistenceProvided <- FALSE;
                xregPersistenceEstimate <- TRUE;
            }
            else if(length(xregPersistence)==1){
                xregPersistence <- rep(xregPersistence,xregNumber);
                xregPersistenceProvided <- TRUE;
                xregPersistenceEstimate <- FALSE;
            }
            else{
                xregPersistenceProvided <- TRUE;
                xregPersistenceEstimate <- FALSE;
            }
        }
        else{
            xregPersistence <- rep(0.05,xregNumber);
            xregPersistenceProvided <- FALSE;
            xregPersistenceEstimate <- TRUE;
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
        xregInitialsProvided <- FALSE;
        xregInitialsEstimate <- FALSE;
        xregPersistenceProvided <- FALSE;
        xregPersistenceEstimate <- FALSE;
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

    #### Checks for the potential number of degrees of freedom ####
    # This is needed in order to make the function work on small samples
    # scale parameter, smoothing parameters and phi
    nParamMax <- (1 + componentsNumber*persistenceEstimate + phiEstimate +
                      # Number of ETS initials
                      (sum(lagsModelAll)-xregNumber-initialNumberARIMA)*(initialType=="optimal") +
                      # ARIMA components: initials + parameters
                      arimaModel*(initialNumberARIMA*(initialType=="optimal") + sum(arOrders) + sum(maOrders)) +
                      # Xreg initials and smoothing parameters
                      xregNumber*(xregInitialsEstimate+xregPersistenceEstimate));

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
            initialNumberARIMA <- componentsNumberARIMA <- 0;
            lagsModelAll <- lagsModelAll[-c(componentsNumber+c(1:componentsNumberARIMA)),,drop=FALSE];
            lagsModelMax <- max(lagsModelAll);

            nParamMax[] <- (1 + componentsNumber*persistenceEstimate + phiEstimate +
                                # Number of ETS initials
                                (sum(lagsModelAll)-xregNumber)*(initialType=="optimal") +
                                # Xreg initials and smoothing parameters
                                xregNumber*(xregInitialsEstimate+xregPersistenceEstimate));
        }
    }

    # If the sample is smaller than the number of parameters
    if(obsNonzero <= nParamMax){
        nParamExo <- xregNumber*(xregInitialsEstimate+xregPersistenceEstimate);
        if(!silent){
            message(paste0("Number of non-zero observations is ",obsNonzero,
                           ", while the maximum number of parameters to estimate is ", nParamMax,".\n",
                           "Updating pool of models."));
        }

        # If the number of observations is still enough for the model selection and the pool is not specified
        if(obsNonzero > (3 + nParamExo) && is.null(modelsPool)){
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
        # If the pool is provided, amend it
        else if(obsNonzero > (3 + nParamExo) & !is.null(modelsPool)){
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
            persistenceEstimate <- FALSE;
            warning("We did not have enough of non-zero observations, so persistence value was set to zero.",
                    call.=FALSE);
            phiEstimate <- FALSE;
        }
        # Can it be even smaller?
        else if(obsNonzero==2){
            modelsPool <- NULL;
            persistence <- 0;
            names(persistence) <- "level";
            persistenceEstimate <- FALSE;
            initialValue <- mean(yInSample);
            initialType <- "provided";
            initialEstimate <- FALSE;
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
            persistenceEstimate <- FALSE;
            initialValue <- yInSample[yInSample!=0];
            initialType <- "p";
            initialEstimate <- FALSE;
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
            persistenceEstimate <- FALSE;
            initialValue <- 0;
            initialType <- "p";
            initialEstimate <- FALSE;
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

    #### Process ellipsis ####
    # Parameters for the optimiser
    if(is.null(ellipsis$maxeval)){
        if(lagsModelMax>12){
            maxeval <- 1000;
        }
        else{
            maxeval <- 200;
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
    if(!any(persistenceEstimate,phiEstimate,initialEstimate,
            xregInitialsEstimate,xregPersistenceEstimate,lambdaEstimate)){
        modelDo <- "use";
    }

    #### Return the values to the previous environment ####
    # Actuals
    assign("y",y,ParentEnvironment);
    assign("yHoldout",yHoldout,ParentEnvironment);
    assign("yInSample",yInSample,ParentEnvironment);
    assign("yNAValues",yNAValues,ParentEnvironment);
    # Index and all related structure variables
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

    # Number of observations and parameters
    assign("obsInSample",obsInSample,ParentEnvironment);
    assign("obsAll",obsAll,ParentEnvironment);
    assign("obsStates",obsStates,ParentEnvironment);
    assign("obsNonzero",obsNonzero,ParentEnvironment);
    assign("obsZero",obsZero,ParentEnvironment);
    assign("parametersNumber",parametersNumber,ParentEnvironment);

    # Model type and lags
    assign("model",model,ParentEnvironment);
    assign("Etype",Etype,ParentEnvironment);
    assign("Ttype",Ttype,ParentEnvironment);
    assign("Stype",Stype,ParentEnvironment);
    assign("modelsPool",modelsPool,ParentEnvironment);
    assign("damped",damped,ParentEnvironment);
    assign("modelDo",modelDo,ParentEnvironment);
    assign("modelIsSeasonal",modelIsSeasonal,ParentEnvironment);
    assign("allowMultiplicative",allowMultiplicative,ParentEnvironment);
    assign("componentsNames",componentsNames,ParentEnvironment);
    assign("componentsNumber",componentsNumber,ParentEnvironment);
    assign("componentsNumberSeasonal",componentsNumberSeasonal,ParentEnvironment);
    # This is the original vector of lags
    assign("lags",lags,ParentEnvironment);
    # This is the vector of lags of ETS components
    assign("lagsModel",lagsModel,ParentEnvironment);
    # This is the vector of all the lags of model (ETS + ARIMA + X)
    assign("lagsModelAll",lagsModelAll,ParentEnvironment);
    # This is the maximum lag
    assign("lagsModelMax",lagsModelMax,ParentEnvironment);

    # Persistence and initials
    assign("persistence",persistence,ParentEnvironment);
    assign("persistenceEstimate",persistenceEstimate,ParentEnvironment);
    assign("phi",phi,ParentEnvironment);
    assign("phiEstimate",phiEstimate,ParentEnvironment);
    assign("initial",initial,ParentEnvironment);
    assign("initialType",initialType,ParentEnvironment);
    assign("initialEstimate",initialEstimate,ParentEnvironment);
    assign("initialValue",initialValue,ParentEnvironment);

    # Occurrence model
    assign("oesModel",oesModel,ParentEnvironment);
    assign("occurrenceModel",occurrenceModel,ParentEnvironment);
    assign("occurrenceModelProvided",occurrenceModelProvided,ParentEnvironment);
    assign("occurrence",occurrence,ParentEnvironment);
    assign("pFitted",pFitted,ParentEnvironment);
    assign("pForecast",pForecast,ParentEnvironment);
    assign("ot",ot,ParentEnvironment);
    assign("otLogical",otLogical,ParentEnvironment);

    # Distribution, loss, bounds and IC
    assign("distribution",distribution,ParentEnvironment);
    assign("loss",loss,ParentEnvironment);
    assign("multisteps",multisteps,ParentEnvironment);
    assign("ic",ic,ParentEnvironment);
    assign("ICFunction",ICFunction,ParentEnvironment);
    assign("bounds",bounds,ParentEnvironment);

    # ARIMA components
    assign("arimaModel",arimaModel,ParentEnvironment);
    assign("arRequired",arRequired,ParentEnvironment);
    assign("iRequired",iRequired,ParentEnvironment);
    assign("maRequired",maRequired,ParentEnvironment);
    assign("nonZeroARI",nonZeroARI,ParentEnvironment);
    assign("nonZeroMA",nonZeroMA,ParentEnvironment);
    assign("lagsModelARIMA",lagsModelARIMA,ParentEnvironment);
    assign("componentsNumberARIMA",componentsNumberARIMA,ParentEnvironment);
    assign("componentsNamesARIMA",componentsNamesARIMA,ParentEnvironment);
    assign("initialNumberARIMA",initialNumberARIMA,ParentEnvironment);
    assign("arOrders",arOrders,ParentEnvironment);
    assign("iOrders",iOrders,ParentEnvironment);
    assign("maOrders",maOrders,ParentEnvironment);

    # Explanatory variables
    assign("xregDo",xregDo,ParentEnvironment);
    assign("xregExist",xregExist,ParentEnvironment);
    assign("xregModel",xregModel,ParentEnvironment);
    assign("xregData",xregData,ParentEnvironment);
    assign("xregNumber",xregNumber,ParentEnvironment);
    assign("xregNames",xregNames,ParentEnvironment);
    assign("xregInitialsProvided",xregInitialsProvided,ParentEnvironment);
    assign("xregInitialsEstimate",xregInitialsEstimate,ParentEnvironment);
    assign("xregPersistenceProvided",xregPersistenceProvided,ParentEnvironment);
    assign("xregPersistenceEstimate",xregPersistenceEstimate,ParentEnvironment);
    assign("xregPersistence",xregPersistence,ParentEnvironment);
    assign("formula",formulaProvided,ParentEnvironment);

    # Ellipsis thingies
    assign("maxeval",maxeval,ParentEnvironment);
    assign("maxtime",maxtime,ParentEnvironment);
    assign("xtol_rel",xtol_rel,ParentEnvironment);
    assign("xtol_abs",xtol_abs,ParentEnvironment);
    assign("algorithm",algorithm,ParentEnvironment);
    assign("print_level",print_level,ParentEnvironment);
    assign("B",B,ParentEnvironment);
    assign("lb",lb,ParentEnvironment);
    assign("ub",ub,ParentEnvironment);
    assign("lambda",lambda,ParentEnvironment);
    assign("lambdaEstimate",lambdaEstimate,ParentEnvironment);
    assign("FI",FI,ParentEnvironment);
}
