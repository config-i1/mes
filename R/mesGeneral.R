parametersChecker <- function(y, model, lags, persistence, phi, initial,
                              distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                                             "dlnorm","dinvgauss"),
                              loss, h, holdout,occurrence,
                              ic=c("AICc","AIC","BIC","BICc"), bounds=c("traditional","admissible","none"),
                              xreg, xregDo, xregInitial, xregPersistence, silent, modelDo, ParentEnvironment,
                              ellipsis, fast=FALSE){

    # The function checks the provided parameters of mes and/or oes
    ##### data #####
    if(any(is.mes.sim(y))){
        y <- y$data;
    }
    else if(inherits(y,"Mdata")){
        h <- y$h;
        holdout <- TRUE;
        if(modelDo!="use"){
            lags <- frequency(y$x);
        }
        y <- ts(c(y$x,y$xx),start=start(y$x),frequency=frequency(y$x));
    }

    if(!is.numeric(y)){
        stop("The provided data is not numeric! Can't construct any model!", call.=FALSE);
    }
    if(!is.null(ncol(y)) && ncol(y)>1){
        # If we deal with data.table, the syntax is different.
        # We don't want to import from data.table, so just use inherits()
        if(inherits(y,"data.table")){
            xreg <- y[,-1];
            y <- y[[1]];
        }
        else{
            xreg <- y[,-1];
            y <- y[,1];
        }
    }

    ####!!! This is a temporary solution !!!####
    # If this is xts, get rid of it
    if(inherits(y,"xts")){
        y <- as.ts(y);
    }

    # Substitute NAs with zeroes.
    ####!!! This will be changed after the introduction of missing data !!!####
    if(any(is.na(y))){
        warning("Data contains NAs. These observations will be substituted by zeroes.",call.=FALSE);
        y[is.na(y)] <- 0;
    }

    # Define obs, the number of observations of in-sample
    obsAll <- length(y) + (1 - holdout)*h;
    obsInSample <- length(y) - holdout*h;
    # dataFreq <- frequency(y);
    # dataStart <- start(y);
    # yForecastStart <- time(y)[obsInSample]+deltat(y);
    yInSample <- matrix(y[1:obsInSample],ncol=1);
    if(holdout){
        yHoldout <- y[-c(1:obsInSample)];
    }
    else{
        yHoldout <- NULL;
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
            # lagsModel <- matrix(lags[1:componentsNumber],ncol=1);
            # lagsModelMax <- max(lagsModel);
            componentsNames <- c(componentsNames[-length(componentsNames)],paste0("seasonal",c(1:(lagsLength-componentsNumber-1))));
            componentsNumberSeasonal[] <- lagsLength-componentsNumber+1;
            # lagsLength <- length(lagsModel);
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
    if(is.oes(occurrence)){
        oesModel <- occurrence;
        occurrence <- oesModel$occurrence;
        occurrenceModelProvided <- TRUE;
    }
    # else if(is.list(occurrence)){
    #     warning(paste0("occurrence is not of the class oes. ",
    #                    "We will try to extract the type of model, but cannot promise anything."),
    #             call.=FALSE);
    #     oesModel <- modelType(occurrence);
    #     occurrence <- occurrence$occurrence;
    #     occurrenceModelProvided <- FALSE;
    # }
    else{
        occurrenceModelProvided <- FALSE;
        oesModel <- NULL;
    }
    pFitted <- matrix(1, obsInSample, 1);

    if(is.numeric(occurrence)){
        # If it is data, then it should correspond to the in-sample.
        if(any(occurrence!=1) && (length(occurrence)!=obsInSample)){
            warning(paste0("Length of the occurrences variable is ",length(occurrence),
                           " when it should be ",obsInSample,".\n",
                           "Switching to occurrence='fixed'."),call.=FALSE);
            occurrence <- "fixed";
        }
        else if(all(occurrence==1)){
            occurrence <- "none";
            occurrenceModelProvided <- FALSE;
        }
        else{
            if(any(occurrence<0,occurrence>1)){
                warning(paste0("Parameter 'occurrence' should contain values between zero and one.\n",
                               "Converting to appropriate vector."),call.=FALSE);
                occurrence[] <- (occurrence!=0)*1;
            }

            # "p" stand for "provided", meaning that we have been provided the values of p
            occurrence <- "provided";
            pFitted[] <- occurrence;
            occurrenceModelProvided <- FALSE;
        }
    }

    occurrence <- match.arg(occurrence[1],c("none","auto","fixed","general","odds-ratio",
                                            "inverse-odds-ratio","direct","provided"));

    otLogical <- (yInSample!=0);

    # If the data is not occurrence, let's assume that the parameter was switched unintentionally.
    if(all(otLogical) & all(occurrence!=c("none","provided"))){
        occurrence <- "none";
        occurrenceModelProvided <- FALSE;
    }

    # This variable just flags, whether we have the occurence in the model or not
    if(occurrence=="none"){
        occurrenceModel <- FALSE;
        otLogical <- rep(TRUE,obsInSample);
    }
    else{
        occurrenceModel <- TRUE;
    }

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

    ot <- matrix(otLogical*1,ncol=1);
    obsNonzero <- sum(ot);
    obsZero <- obsInSample - obsNonzero;

    # Update the number of parameters
    if(occurrenceModelProvided){
        parametersNumber[2,3] <- nparam(oesModel);
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

    #### Explanatory variables: xreg, xregDo, xregInitial, xregPersistence ####
    xregDo <- match.arg(xregDo,c("use","select"));
    xregExist <- !is.null(xreg);
    if(!xregExist){
        xregDo <- "use";
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
            xregInitialiser <- function(Etype,distribution){
                # Fix the default distribution for ALM
                if(distribution=="default"){
                    distribution <- switch(Etype,
                                           "A"="dnorm",
                                           "M"="dlnorm");
                }
                # Return the estimated model based on the provided xreg
                if(Etype=="M" && any(distribution==c("dnorm","dlogis","dlaplace","dt","ds","dalaplace"))){
                    return(alm(log(y)~.,xregData,distribution=distribution));
                }
                else{
                    return(alm(y~.,xregData,distribution=distribution));
                }
            }
            # Extract names and form a proper matrix for the regression
            xregNames <- c("y",colnames(xreg));
            xregData <- cbind(yInSample,xreg[1:obsInSample,]);
            colnames(xregData) <- xregNames;

            if(Etype!="Z"){
                testModel <- xregInitialiser(Etype,distribution);
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
                testModel <- xregInitialiser("A",distribution);
                xregModel[[1]]$xregInitial <- testModel$coefficients[-1];
                xregModel[[1]]$other <- testModel$other;
                # Multiplicative model
                testModel[] <- xregInitialiser("M",distribution);
                xregModel[[2]]$xregInitial <- testModel$coefficients[-1];
                xregModel[[2]]$other <- testModel$other;
            }

            # Write down the number and names of parameters
            xregNumber <- ncol(testModel$data)-1;
            xregNames <- colnames(testModel$data)[-1];
            xregData <- testModel$data[,-1];
            if(nrow(xreg)>obsInSample){
                xregData <- as.matrix(model.frame(~.,data=xreg))[,xregNames];
            }
        }
        else{
            xregInitialsProvided <- TRUE;
            xregInitialsEstimate <- FALSE;

            xregModel[[1]]$xregInitial <- xregInitial;
            if(Etype=="Z"){
                xregModel[[2]]$xregInitial <- xregInitial;
            }

            # Write down the number and names of parameters
            if(nrow(xreg>obsAll)){
                xregData <- xreg[1:obsAll,];
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
        lagsModelAll <- matrix(c(lagsModel,rep(1,xregNumber)),ncol=1);
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
        lagsModelAll <- lagsModel;
    }

    #### Process ellipsis ####
    # Parameters for the optimiser
    if(is.null(ellipsis$maxeval)){
        maxeval <- 500;
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
    assign("lags",lags,ParentEnvironment);
    assign("lagsModel",lagsModel,ParentEnvironment);
    assign("lagsModelMax",lagsModelMax,ParentEnvironment);
    assign("lagsModelAll",lagsModelAll,ParentEnvironment);
    assign("lagsLength",lagsLength,ParentEnvironment);

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
    assign("ot",ot,ParentEnvironment);
    assign("otLogical",otLogical,ParentEnvironment);

    # Distribution, loss, bounds and IC
    assign("distribution",distribution,ParentEnvironment);
    assign("loss",loss,ParentEnvironment);
    assign("multisteps",multisteps,ParentEnvironment);
    assign("ic",ic,ParentEnvironment);
    assign("ICFunction",ICFunction,ParentEnvironment);
    assign("bounds",bounds,ParentEnvironment);

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
