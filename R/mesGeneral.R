parametersChecker <- function(y, model, lags, date, persistence, phi, initia, loss, distribution, occurrence, ic, bounds,
                              xreg, xregDo, xregInitial, xregPersistence, silent, fast, ParentEnvironment, ...){

    # The function checks the provided parameters of mes and/or omes
    ##### data #####
    if(any(is.mes.sim(y))){
        y <- y$data;
    }
    else if(inherits(y,"Mdata")){
        # h <- y$h;
        # holdout <- TRUE;
        # y <- ts(c(y$x,y$xx),start=start(y$x),frequency=frequency(y$x));
        y <- y$x;
    }

    if(!is.numeric(y)){
        stop("The provided data is not numeric! Can't construct any model!", call.=FALSE);
    }
    if(!is.null(ncol(y))){
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

    # Substitute NAs with zeroes.
    ####!!! This will be changed after the introduction of missing data !!!####
    if(any(is.na(y))){
        warning("Data contains NAs. These observations will be substituted by zeroes.",call.=FALSE);
        y[is.na(y)] <- 0;
    }

    # Define obs, the number of observations of in-sample
    obsInSample <- length(y);
    # dataFreq <- frequency(y);
    # dataStart <- start(y);
    # yForecastStart <- time(y)[obsInSample]+deltat(y);

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
        else if(any(unlist(strsplit(model,""))=="Z") |
                any(unlist(strsplit(model,""))=="X") |
                any(unlist(strsplit(model,""))=="Y")){
            modelDo <- "select";
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
        warning(paste0("Wrong error type: ",Etype,". Should be 'Z', 'X', 'Y', 'A' or 'M'.\n",
                       "Changing to 'Z'"),call.=FALSE);
        Etype <- "Z";
        modelDo <- "select";
    }

    ### Check trend type
    if(all(Ttype!=c("Z","X","Y","N","A","M","C"))){
        warning(paste0("Wrong trend type: ",Ttype,". Should be 'Z', 'X', 'Y', 'N', 'A' or 'M'.\n",
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

    # Get rid of duplicates in lags
    if(length(unique(lags))!=length(lags)){
        lags <- unique(lags);
    }

    # Lags of the model
    lagsModel <- matrix(lags,ncol=1);
    lagsModelMax <- max(lagsModel);
    lagsLength <- length(lags);

    modelIsSeasonal <- Stype!="N";
    #### Check the seasonal model vs lags ####
    if(all(Stype!=c("Z","X","Y","N","A","M","C"))){
        warning(paste0("Wrong seasonality type: ",Stype,". Should be 'Z', 'X', 'Y', 'C', 'N', 'A' or 'M'.",
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
    }

    # Check, whether the number of lags and the number of components are the same
    if(lagsLength>componentsNumber){
        if(Stype!="N"){
            componentsNames <- c(componentsNames[-length(componentsNames)],paste0("seasonal",c(1:(lagsLength-componentsNumber))));
            componentsNumber[] <- lagsLength;
        }
        else{
            lagsModel <- matrix(lags[1:componentsNumber],ncol=1);
            lagsModelMax <- max(lagsModel);
            lagsLength <- length(lags);
        }
    }
    else if(lagsLength<componentsNumber){
        stop("The number of components of the model is smaller than the number of provided lags", call.=FALSE);
    }

    #### Observations in the states matrix ####
    # Define the number of rows that should be in the matvt
    obsStates <- obsInSample + 2*lagsModelMax;

    #### Check the dates ####

    #### Distribution selected ####
    distribution <- match.arg(distribution,c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                                             "dlnorm","dbcnorm","dinvgauss"));

    #### Loss function type ####
    loss <- match.arg(loss,c("likelihood","MSE","MAE","HAM",
                             "MSEh","TMSE","GTMSE","MSCE",
                             "MAEh","TMAE","GTMAE","MACE",
                             "HAMh","THAM","GTHAM","CHAM",
                             "GPL","aMSEh","aTMSE","aGTMSE","aGPL"));

    if(any(loss==c("MSEh","TMSE","GTMSE","MSCE","MAEh","TMAE","GTMAE","MACE",
                     "HAMh","THAM","GTHAM","CHAM",
                     "GPL","aMSEh","aTMSE","aGTMSE","aGPL"))){
        multisteps <- TRUE;
    }
    else{
        multisteps <- FALSE;
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
            phi <- NULL;
            phiEstimate <- TRUE;
        }
        else if(is.numeric(phi) & (phi<0 | phi>2)){
            warning(paste0("Damping parameter should lie in (0, 2) region. ",
                           "Changing to the estimation of phi."),call.=FALSE);
            phi <- NULL;
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
        }
        else{
            phiEstimate <- FALSE;
        }
    }

    #### Vector of initial values ####
    # initial type can be: "o" - optimal, "b" - backcasting, "p" - provided.
    if(is.character(initial)){
        initialType <- match.arg(initial, c("optimal","backcasting"));
    }
    else if(is.null(initial)){
        if(silentText){
            message("Initial value is not selected. Switching to optimal.");
        }
        initialType <- "optimal";
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
                    warning(paste0("Wrong length of initial vector. Should be ",sum(lags),
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

    #### Occurrence variable ####

    #### Information Criteria ####
    ic <- ic[1];
    if(all(ic!=c("AICc","AIC","BIC","BICc"))){
        warning(paste0("Strange type of information criteria defined: ",ic,". Switching to 'AICc'."),
                call.=FALSE);
        ic <- "AICc";
    }

    #### Bounds for the smoothing parameters ####
    bounds <- match.arg(bounds,c("usual","admissible","none"));

    #### Explanatory variables: xreg, xregDo, xregInitial, xregPersistence ####

}
