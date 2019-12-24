parametersChecker <- function(...){
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

}
