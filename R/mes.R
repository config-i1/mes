#' Mixed Exponential Smoothing
#'
#' Function constructs an advanced Single Source of Error model, based on ETS
#' taxonomy.
#'
#' Function estimates ETS in a form of the Single Source of Error state space
#' model of the following type:
#'
#' \deqn{y_{t} = o_t (w(v_{t-l}) + h(x_t, a_{t-1}) + r(v_{t-l}) \epsilon_{t})}
#'
#' \deqn{v_{t} = f(v_{t-l}, a_{t-1}) + g(v_{t-l}, a_{t-1}, x_{t}) \epsilon_{t}}
#'
#' Where \eqn{o_{t}} is the Bernoulli distributed random variable (in case of
#' normal data it equals to 1 for all observations), \eqn{v_{t}} is the state
#' vector and \eqn{l} is the vector of lags, \eqn{x_t} is the vector of
#' exogenous variables. w(.) is the measurement function, r(.) is the error
#' function, f(.) is the transition function, g(.) is the persistence
#' function and \eqn{a_t} is the vector of parameters for exogenous variables.
#' Finally, \eqn{\epsilon_{t}} is the error term.
#'
#' The implemented model allows introducing several seasonal states and supports
#' intermittent data via the \code{occurrence} variable.
#'
#' The error term \eqn{\epsilon_t} can follow different distributions, which
#' are regulated via the \code{distribution} parameter. This includes:
#' \enumerate{
#' \item \code{default} - Normal distribution is used for the Additive error models,
#' Inverse Gaussian is used for the Multiplicative error models.
#' \item \link[stats]{dnorm} - Normal distribution,
#' \item \link[stats]{dlogis} - Logistic Distribution,
#' \item \link[greybox]{dlaplace} - Laplace distribution,
#' \item \link[greybox]{dalaplace} - Asymmetric Laplace distribution,
#' \item \link[stats]{dt} - T-distribution,
#' \item \link[greybox]{ds} - S-distribution,
#' \item \link[stats]{dlnorm} - Log normal distribution,
# \item \link[greybox]{dbcnorm} - Box-Cox normal distribution,
#' \item \link[statmod]{dinvgauss} - Inverse Gaussian distribution,
#' }
#'
#' For some more information about the model and its implementation, see the
#' vignette: \code{vignette("mes","mes")}.
#'
#' @template ssAuthor
#' @template ssKeywords
#'
#' @template ssGeneralRef
#' @template ssIntermittentRef
#' @template ssETSRef
#' @template ssIntervalsRef
#'
#' @param y Vector, containing data needed to be forecasted. If a matrix is
#' provided, then the first column is used as a response variable, while the rest
#' of the matrix is used as a set of explanatory variables.
#' @param model The type of ETS model. The first letter stands for the type of
#' the error term ("A" or "M"), the second (and sometimes the third as well) is for
#' the trend ("N", "A", "Ad", "M" or "Md"), and the last one is for the type of
#' seasonality ("N", "A" or "M"). In case of several lags, the seasonal components
#' are assumed to be the same. The model is then printed out as
#' MES(M,Ad,M[m1],M[m2],...), where m1, m2, ... are the lags specified by the
#' \code{lags} parameter.
#'
#' \code{ZZZ} means that the model will be selected based on the
#' chosen information criteria type. Models pool can be restricted with additive
#' only components via \code{model="XXX"}. Selection between multiplicative models
#' (excluding additive components) is regulated using \code{model="YYY"}. Finally,
#' \code{model="CCC"} trigers the combination of forecasts of models using AIC
#' weights (Kolassa, 2011). All of this can be finely tuned. For example,
#' \code{model="CCN"} will combine forecasts of all non-seasonal models and
#' \code{model="CXY"} will combine forecasts of all the models with
#' non-multiplicative trend and non-additive seasonality with either additive
#' or multiplicative error.
#'
#' The parameter \code{model} can also be a vector of names of models for a
#' finer tuning (pool of models). For example, \code{model=c("ANN","AAA")} will
#' estimate only two models and select the best of them.
#'
#' Also \code{model} can accept a previously estimated mes model and use all
#' its parameters.
#'
#' Keep in mind that model selection with "Z" components uses Branch and Bound
#' algorithm and may skip some models that could have slightly smaller
#' information criteria. If you want to do a exhaustive search, you would need
#' to list all the models to check as a vector.
#' @param lags Defines lags for the corresponding components. All components
#' count, starting from level, so ETS(M,M,M) model for monthly data will have
#' lags=c(1,1,12). If fractional numbers are provided, then it is assumed that
#' the data is not periodic. The parameter \code{date} is then needed in order
#' to setup the appropriate time series structure.
#' @param distribution what density function to assume for the error term. The full
#' name of the distribution should be provided, starting with the letter "d" -
#' "density". The names align with the names of distribution functions in R.
#' For example, see \link[stats]{dnorm}.
#' @param loss The type of Loss Function used in optimization. \code{loss} can
#' be:
#' \itemize{
#' \item \code{likelihood} - the model is estimated via the maximisation of the
#' likelihood of the function specified in \code{distribution};
#' \item \code{MSE} (Mean Squared Error),
#' \item \code{MAE} (Mean Absolute Error),
#' \item \code{HAM} (Half Absolute Moment),
#' \item \code{LASSO} - use LASSO to shrink the parameters of the model;
#' \item \code{RIDGE} - use RIDGE to shrink the parameters of the model;
#' \item \code{TMSE} - Trace Mean Squared Error,
#' \item \code{GTMSE} - Geometric Trace Mean Squared Error,
#' \item \code{MSEh} - optimisation using only h-steps ahead error,
#' \item \code{MSCE} - Mean Squared Cumulative Error.
#' }
#' Note that model selection and combination works properly only for the default
#' \code{loss="likelihood"}.
#'
#' There are also available analytical approximations for multistep functions:
#' \code{aMSEh}, \code{aTMSE} and \code{aGTMSE}. These can be useful in cases
#' of small samples.
#'
#' Finally, just for fun the absolute and half analogues of multistep estimators
#' are available: \code{MAEh}, \code{TMAE}, \code{GTMAE}, \code{MACE}, \code{TMAE},
#' \code{HAMh}, \code{THAM}, \code{GTHAM}, \code{CHAM}.
#' @param h The forecast horizon. Mainly needed for the multistep loss functions.
#' @param holdout Logical. If \code{TRUE}, then the holdout of the size \code{h}
#' is taken from the data (can be used for the model testing purposes.
#' @param persistence Persistence vector \eqn{g}, containing smoothing
#' parameters. If \code{NULL}, then estimated.
#' @param phi Value of damping parameter. If \code{NULL} then it is estimated.
#' Only applicable for damped-trend models.
#' @param initial Can be either character or a vector of initial states. If it
#' is character, then it can be \code{"optimal"}, meaning that the initial
#' states are optimised, or \code{"backcasting"}, meaning that the initials are
#' produced using backcasting procedure (advised for data with high frequency).
#' @param occurrence The type of model used in probability estimation. Can be
#' \code{"none"} - none,
#' \code{"fixed"} - constant probability,
#' \code{"general"} - the general Beta model with two parameters,
#' \code{"odds-ratio"} - the Odds-ratio model with b=1 in Beta distribution,
#' \code{"inverse-odds-ratio"} - the model with a=1 in Beta distribution,
#' \code{"direct"} - the TSB-like (Teunter et al., 2011) probability update
#' mechanism a+b=1,
#' \code{"auto"} - the automatically selected type of occurrence model.
#'
#' The type of model used in the occurrence is equal to the one provided in the
#' \code{model} parameter.
#'
#' Also, a model produced using omess() (NOT AVAILABLE YET), \link[smooth]{oes} or
#' \link[greybox]{alm} function can be used here.
#' @param ic The information criterion to use in the model selection / combination
#' procedure.
#' @param bounds The type of bounds for the persistence to use in the model
#' estimation. Can be either \code{admissible} - guaranteeing the stability of the
#' model, \code{traditional} - restricting the values with (0, 1) or \code{none} - no
#' restrictions (potentially dangerous).
#' @param silent Specifies, whether to provide the progress of the function or not.
#' If \code{TRUE}, then the function will print what it does and how much it has
#' already done.
#' @param xreg The vector (either numeric or time series) or the matrix (or
#' data.frame / data.table) of exogenous variables that should be included in the
#' model. If matrix is included than columns should contain variables and rows -
#' observations.
#' Note that \code{xreg} should have number of observations equal to
#' the length of the response variable \code{y}. If it is not equal, then the
#' function will either trim or extrapolate the data.
#' @param xregDo The variable defines what to do with the provided xreg:
#' \code{"use"} means that all of the data should be used, while
#' \code{"select"} means that a selection using \code{ic} should be done.
#' @param xregInitial The vector of initial parameters for exogenous variables.
#' Ignored if \code{xreg} is NULL.
#' @param xregPersistence The persistence vector \eqn{g_X}, containing smoothing
#' parameters for exogenous variables. If \code{NULL}, then estimated. If \code{0}
#' then each element of the vector is set to zero. Prerequisite - non-NULL \code{xreg}.
# @param fast if \code{TRUE}, then the function won't check whether
#' the provided vectors are correct and will use them directly in the model
#' construction.
#' @param ...  Other non-documented parameters. For example \code{FI=TRUE} will
#' make the function also produce Fisher Information matrix, which then can be
#' used to calculated variances of smoothing parameters and initial states of
#' the model. This is used in the \link[stats]{vcov} method.
#' Starting values of parameters can be passed via \code{B}, while the upper and lower
#' bounds should be passed in \code{ub} and \code{lb} respectively. In this case they
#' will be used for optimisation. These values should have the length equal
#' to the number of parameters to estimate in the following order:
#' \enumerate{
#' \item All smoothing parameters (for the states and then for the explanatory variables);
#' \item Damping parameter (if needed);
#' \item All the initial values (for the states and then for the explanatory variables).
#' }
#' You can also pass parameters to the optimiser in order to fine tune its work:
#' \itemize{
#' \item \code{maxeval} - maximum number of evaluations to carry out (default is 100);
#' \item \code{maxtime} - stop, when the optimisation time (in seconds) exceeds this;
#' \item \code{xtol_rel} - the precision of the optimiser (the default is 1E-6);
#' \item \code{algorithm} - the algorithm to use in optimisation
#' (\code{"NLOPT_LN_BOBYQA"} by default);
#' \item \code{print_level} - the level of output for the optimiser (0 by default).
#' }
#' You can read more about these parameters by running the function
#' \link[nloptr]{nloptr.print.options()}.
#' Finally, the parameter \code{lambda} for LASSO / RIDGE, Asymmetric Laplace and df
#' of Student's distribution can be provided here as well.
#'
#' @return Object of class "mes" is returned. It contains the list of the
#' following values:
#' \itemize{
#' \item \code{model} - the name of the constructed model,
#' \item \code{timeElapsed} - the time elapsed for the estimation of the model,
#' \item \code{y} - the in-sample part of the data used for the training of the model,
#' \item \code{holdout} - the holdout part of the data, excluded for purposes of model evaluation,
#' \item \code{fitted} - the vector of fitted values,
#' \item \code{residuals} - the vector of residuals,
#' \item \code{forecast} - the point forecast for h steps ahead (by default NA is returned),
#' \item \code{states} - the matrix of states with observations in rows and states in columns,
#' \item \code{persisten} - the vector of smoothing parameters,
#' \item \code{phi} - the value of damping parameter,
#' \item \code{transition} - the transition matrix,
#' \item \code{measurement} - the measurement matrix with observations in rows and state elements
#' in columns,
#' \item \code{initialType} - the type of initialisation used ("optimal" / "backcasting" / "provided"),
#' \item \code{initial} - the initial values, including level, trend and seasonal components,
#' \item \code{nParam} - the matrix of the estimated / provided parameters,
#' \item \code{occurrence} - the oes model used for the occurrence part of the model,
#' \item \code{xreg} - the matrix of explanatory variables after all expansions and transformations,
#' \item \code{xregInitial} - the vector of initials for the parameters of explanatory variable,
#' \item \code{xregPersistence} - the vector of smoothing parameters for the explanatory variables,
#' \item \code{loss} - the type of loss function used in the estimation,
#' \item \code{lossValue} - the value of that loss function,
#' \item \code{logLik} - the value of the log-likelihood,
#' \item \code{distribution} - the distribution function used in the calculation of the likelihood,
#' \item \code{scale} - the value of the scale parameter,
#' \item \code{lambda} - the value of the parameter used in LASSO / dalaplace / dt,
#' \item \code{B} - the vector of all estimated parameters,
#' \item \code{lags} - the vector of lags used in the model construction.
#' }
#'
#' @seealso \code{\link[forecast]{ets}, \link[smooth]{es}}
#'
#' @examples
#'
#' # Model selection using a specified pool of models
#' ourModel <- mes(rnorm(100,100,10), model=c("ANN","ANA","AAA"), lags=c(5,10))
#'
#' summary(ourModel)
#' forecast(ourModel)
#' plot(forecast(ourModel))
#'
#' # Model combination using a specified pool
#' ourModel <- mes(rnorm(100,100,10), model=c("ANN","AAN","MNN","CCC"), lags=c(5,10))
#'
#' @importFrom forecast forecast
#' @importFrom greybox dlaplace dalaplace ds stepwise alm
#' @importFrom stats dnorm dlogis dt dlnorm frequency
#' @importFrom statmod dinvgauss
#' @importFrom nloptr nloptr
#' @importFrom pracma hessian
#' @useDynLib mes
#' @export mes
mes <- function(y, model="ZZZ", lags=c(frequency(y)),
                distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                               "dlnorm","dinvgauss"),
                loss=c("likelihood","MSE","MAE","HAM","LASSO","RIDGE","MSEh","TMSE","GTMSE","MSCE"),
                h=0, holdout=FALSE,
                persistence=NULL, phi=NULL, initial=c("optimal","backcasting"),
                occurrence=c("none","auto","fixed","general","odds-ratio","inverse-odds-ratio","direct"),
                ic=c("AICc","AIC","BIC","BICc"), bounds=c("usual","admissible","none"),
                xreg=NULL, xregDo=c("use","select"), xregInitial=NULL, xregPersistence=0,
                silent=TRUE, ...){
    # Copyright (C) 2019 - Inf  Ivan Svetunkov
    #
    # Parameters that were moved to forecast() and predict() functions:
    # h=10, holdout=FALSE, cumulative=FALSE,
    # interval=c("none","parametric","likelihood","semiparametric","nonparametric","confidence"), level=0.95,

    # Start measuring the time of calculations
    startTime <- Sys.time();

    ellipsis <- list(...);
    # If a previous model is provided as a model, write down the variables
    if(is.mes(model) || is.mes.sim(model)){
        # If this is the simulated data, extract the parameters
        # if(is.mes.sim(model) & !is.null(dim(model$data))){
        #     warning("The provided model has several submodels. Choosing a random one.",call.=FALSE);
        #     i <- round(runif(1,1:length(model$persistence)));
        #     persistence <- model$persistence[,i];
        #     initial <- model$initial[,i];
        #     initialSeason <- model$initialSeason[,i];
        #     if(any(model$iprob!=1)){
        #         occurrence <- "a";
        #     }
        # }
        # else{
        persistence <- model$persistence;
        initial <- model$initial;
        initialSeason <- model$initialSeason;
        #     if(any(model$iprob!=1)){
        #         occurrence <- "a";
        #     }
        # }
        lags <- model$lags;
        distribution <- model$distribution;
        loss <- model$loss;
        persistence <- model$persistence;
        phi <- model$phi;
        if(model$initialType!="backcasting"){
            initial <- model$initial;
        }
        else{
            initial <- "b";
        }
        occurrence <- model$occurrence;
        ic <- model$ic;
        bounds <- model$bounds;
        lambda <- model$lambda;
        ellipsis$B <- model$B;
        CFValue <- model$lossValue;
        logLikMESValue <- logLik(model);
        if(is.null(xreg)){
            xreg <- model$xreg;
        }
        else{
            if(is.null(model$xreg)){
                xreg <- NULL;
            }
            else{
                if(ncol(xreg)!=ncol(model$xreg)){
                    xreg <- xreg[,colnames(model$xreg)];
                }
            }
        }

        xregInitial <- model$xregInitial;
        xregPersistence <- model$xregPersistence;

        model <- modelType(model);
        modelDo <- "use";
        # if(any(unlist(gregexpr("C",model))!=-1)){
        #     initial <- "o";
        # }
    }
    else if(inherits(model,"ets")){
        # Extract smoothing parameters
        i <- 1;
        lags <- 1;
        persistence <- coef(model)[i];
        if(model$components[2]!="N"){
            i <- i+1;
            persistence <- c(persistence,coef(model)[i]);
            if(model$components[3]!="N"){
                i <- i+1;
                persistence <- c(persistence,coef(model)[i]);
            }
        }
        else{
            if(model$components[3]!="N"){
                i <- i+1;
                persistence <- c(persistence,coef(model)[i]);
            }
        }

        # Damping parameter
        if(model$components[4]=="TRUE"){
            i <- i+1;
            phi <- coef(model)[i];
        }

        # Initials
        i <- i+1;
        initial <- coef(model)[i];
        # Initial for the trend
        if(model$components[2]!="N"){
            i <- i+1;
            lags <- c(lags,1);
            initial <- c(initial,coef(model)[i]);
        }

        # Initials of seasonal component
        if(model$components[3]!="N"){
            if(model$components[2]!="N"){
                initial <- c(initial,rev(model$states[1,-c(1:2)]));
            }
            else{
                initial <- c(initial,rev(model$states[1,-c(1)]));
            }
            lags <- c(lags,model$m);
        }
        model <- modelType(model);
        distribution <- "dnorm";
        loss <- "likelihood";
        modelDo <- "use"
    }
    else if(is.character(model)){
        modelDo <- "";
        # Everything is okay
    }
    else{
        modelDo <- "";
        warning("A model of an unknown class was provided. Switching to 'ZZZ'.",call.=FALSE);
        model <- "ZZZ";
    }

    #### Check the parameters of the function and create variables based on them ####
    parametersChecker(y, model, lags, persistence, phi, initial,
                      distribution, loss, h, holdout, occurrence, ic, bounds,
                      xreg, xregDo, xregInitial, xregPersistence,
                      silent, modelDo, ParentEnvironment=environment(), ellipsis);

    #### The function creates the technical variables (lags etc) based on the type of the model ####
    architector <- function(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType){
        if(Ttype!="N"){
            # Make lags (1, 1) if they are not
            lagsModel <- matrix(c(1,1),ncol=1);
            componentsNames <- c("level","trend");
        }
        else{
            # Make lags (1, ...) if they are not
            lagsModel <- matrix(c(1),ncol=1);
            componentsNames <- c("level");
        }
        if(Stype!="N"){
            # If the lags are for the non-seasonal model
            lagsModel <- matrix(c(lagsModel,lags[lags>1]),ncol=1);
            componentsNumberSeasonal <- sum(lags>1);
            if(componentsNumberSeasonal>1){
                componentsNames <- c(componentsNames,paste0("seasonal",c(1:componentsNumberSeasonal)));
            }
            else{
                componentsNames <- c(componentsNames,"seasonal");
            }
        }
        else{
            componentsNumberSeasonal <- 0;
        }

        componentsNumber <- lagsLength <- length(lagsModel);
        lagsModelAll <- matrix(c(lagsModel,rep(1,xregNumber)),ncol=1);
        lagsModelMax <- max(lagsModelAll);

        # Define the number of cols that should be in the matvt
        obsStates <- obsInSample + lagsModelMax*switch(initialType,
                                                       "backcasting"=2,
                                                       1);

        return(list(lagsModel=lagsModel,lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax, lagsLength=lagsLength,
                    componentsNumber=componentsNumber, componentsNumberSeasonal=componentsNumberSeasonal,
                    componentsNames=componentsNames, obsStates=obsStates));
    }

    #### The function creates the necessary matrices based on the model and provided parameters ####
    # This is needed in order to initialise the estimation
    creator <- function(Etype, Ttype, Stype,
                        lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                        obsStates, obsInSample, componentsNumber, componentsNumberSeasonal,
                        componentsNames, otLogical,
                        yInSample, persistence, persistenceEstimate, phi,
                        initialValue, initialType,
                        xregProvided, xregInitialsProvided, xregPersistence,
                        xregModel, xregData, xregNumber, xregNames){
        # Matrix of states. Time in columns, components in rows
        matVt <- matrix(NA, componentsNumber+xregNumber, obsStates, dimnames=list(c(componentsNames,xregNames),NULL));

        # Measurement rowvector
        matWt <- matrix(1, obsInSample, componentsNumber+xregNumber, dimnames=list(NULL,c(componentsNames,xregNames)));
        # If xreg are provided, then fill in the respective values in Wt vector
        if(xregProvided){
            matWt[,componentsNumber+1:xregNumber] <- xregData;
        }

        # Transition matrix
        matF <- diag(componentsNumber+xregNumber);

        # Persistence vector
        vecG <- matrix(0, componentsNumber+xregNumber, 1, dimnames=list(c(componentsNames,xregNames),NULL));
        if(!persistenceEstimate){
            vecG[1:componentsNumber,] <- persistence;
        }
        if(xregProvided){
            vecG[componentsNumber+1:xregNumber,] <- xregPersistence;
        }

        # Damping parameter value
        if(Ttype!="N"){
            matF[1,2] <- phi;
            matF[2,2] <- phi;

            matWt[,2] <- phi;
        }

        # Calculate the initials for the matVt and insert them
        if(initialType!="provided"){
            # For the seasonal models
            if(Stype!="N"){
                yDecomposition <- msdecompose(yInSample, lags[lags!=1], type=switch(Stype, "A"="additive", "M"="multiplicative"));
                j <- 1;
                # level
                matVt[j,1:lagsModelMax] <- mean(yInSample[1:lagsModelMax]);
                j <- j+1;
                if(Ttype!="N"){
                    if(Ttype=="A" && Stype=="M"){
                        matVt[j,1:lagsModelMax] <- prod(yDecomposition$initial)-yDecomposition$initial[1];
                    }
                    else if(Ttype=="M" && Stype=="A"){
                        matVt[j,1:lagsModelMax] <- sum(yDecomposition$initial)/yDecomposition$initial[1];
                    }
                    else{
                        matVt[j,1:lagsModelMax] <- yDecomposition$initial[2];
                    }
                    j <- j+1;
                }
                for(i in 1:componentsNumberSeasonal){
                    matVt[i+j-1,(lagsModelMax-lagsModel[i+j-1])+1:lagsModel[i+j-1]] <- yDecomposition$seasonal[[i]];
                }
            }
            # Non-seasonal models
            else{
                matVt[1,1] <- mean(yInSample[1:max(lagsModelMax,obsInSample*0.2)]);
                if(Ttype!="N"){
                    matVt[2,1] <- switch(Ttype,
                                         "A" = mean(diff(yInSample),na.rm=TRUE),
                                         "M" = exp(mean(diff(log(yInSample[otLogical])),na.rm=TRUE)));
                }
            }
            if(Etype=="M" && any(matVt[1,1:lagsModelMax]==0)){
                matVt[1,1:lagsModelMax] <- mean(yInSample);
            }
        }
        # Else, insert the provided ones
        else{
            j <- 1;
            matVt[j,1:lagsModelMax] <- initialValue[j];
            j <- j+1;
            if(Ttype!="N"){
                matVt[j,1:lagsModelMax] <- initialValue[j];
                j <- j+1;
            }
            if(Stype!="N"){
                for(i in j:componentsNumber){
                    indices <- sum(lagsModel[1:(j-1)])+1:lagsModel[i];
                    matVt[i,(lagsModelMax-lagsModel[i])+1:lagsModel[i]] <- initialValue[indices];
                }
            }
        }

        # Fill in the initials for xreg
        if(xregExist){
            if(Etype=="A"){
                matVt[componentsNumber+1:xregNumber,1:lagsModelMax] <- xregModel[[1]]$xregInitial;
            }
            else{
                matVt[componentsNumber+1:xregNumber,1:lagsModelMax] <- xregModel[[2]]$xregInitial;
            }
        }

        return(list(matVt=matVt, matWt=matWt, matF=matF, vecG=vecG));
    }

    #### The function fills in the existing matrices with values of A ####
    # This is needed in order to do the estimation and the fit
    filler <- function(B,
                       Ttype, Stype, componentsNumber, lagsModel, lagsModelMax,
                       matVt, matWt, matF, vecG,
                       persistenceEstimate, phiEstimate, initialType,
                       xregInitialsEstimate, xregPersistenceEstimate, xregNumber){
        j <- 1;
        # Fill in persistence
        if(persistenceEstimate){
            vecG[j:componentsNumber] <- B[j:componentsNumber];
            j <- j+componentsNumber;
        }

        # Persistence if xreg is provided
        if(xregPersistenceEstimate){
            vecG[j+1:xregNumber-1] <- B[j+1:xregNumber-1];
            j <- j+xregNumber;
        }

        # Damping parameter
        if(phiEstimate){
            matWt[,2] <- B[j];
            matF[1:2,2] <- B[j];
            j <- j+1;
        }

        # Initials
        if(initialType=="optimal"){
            i <- 1;
            matVt[i,1:lagsModelMax] <- B[j];
            j <- j+1;
            i <- i+1;
            if(Ttype!="N"){
                matVt[i,1:lagsModelMax] <- B[j];
                j <- j+1;
                i <- i+1;
            }
            if(Stype!="N"){
                for(k in i:componentsNumber){
                    matVt[k,(lagsModelMax-lagsModel[k])+1:lagsModel[k]] <- B[j+0:(lagsModel[k]-1)];
                    j <- j+lagsModel[k];
                }
            }
        }

        # Initials of the xreg
        if(xregInitialsEstimate){
            matVt[componentsNumber+1:xregNumber,1:lagsModelMax] <- B[j+1:xregNumber-1];
        }

        return(list(matVt=matVt, matWt=matWt, matF=matF, vecG=vecG));
    }

    #### The function initialises the vector B for ETS ####
    initialiser <- function(Etype, Ttype, Stype, componentsNumberSeasonal,
                            componentsNumber, lagsModel, lagsModelMax, matVt,
                            persistenceEstimate, phiEstimate, initialType,
                            xregInitialsEstimate, xregPersistenceEstimate, xregNumber, lambdaEstimate){
        # Persistence of states, persistence of xreg, phi, initials, initials for xreg
        B <- Bl <- Bu <- vector("numeric",persistenceEstimate*componentsNumber+xregPersistenceEstimate*xregNumber+phiEstimate+
                                    (initialType=="optimal")*sum(lagsModel)+xregInitialsEstimate*xregNumber+lambdaEstimate);

        j <- 1;
        # Fill in persistence
        if(persistenceEstimate){
            if(any(c(Etype,Ttype,Stype)=="M")){
                B[j:componentsNumber] <- c(0.01,0.005,rep(0.01,componentsNumberSeasonal))[j:componentsNumber];
            }
            else{
                B[j:componentsNumber] <- c(0.1,0.05,rep(0.1,componentsNumberSeasonal))[j:componentsNumber];
            }
            Bl[j:componentsNumber] <- rep(-5, componentsNumber);
            Bu[j:componentsNumber] <- rep(5, componentsNumber);
            names(B)[1] <- "alpha";
            if(Ttype!="N"){
                names(B)[2] <- "beta";
            }
            if(Stype!="N"){
                if(componentsNumberSeasonal>1){
                    names(B)[(componentsNumber-componentsNumberSeasonal+1):componentsNumber] <-
                        paste0("gamma",c(1:componentsNumberSeasonal));
                }
                else{
                    names(B)[(componentsNumber-componentsNumberSeasonal+1):componentsNumber] <- "gamma";
                }
            }
            j <- j+componentsNumber;
        }

        # Persistence if xreg is provided
        if(xregPersistenceEstimate){
            B[j-1+1:xregNumber] <- rep(switch(Etype,
                                              "A"=0.1,
                                              "M"=0.01),xregNumber);
            Bl[j-1+1:xregNumber] <- rep(-5, xregNumber);
            Bu[j-1+1:xregNumber] <- rep(5, xregNumber);
            names(B)[j-1+1:xregNumber] <- paste0("delta",c(1:xregNumber));
            j <- j+xregNumber;
        }

        # Damping parameter
        if(phiEstimate){
            B[j] <- 0.95;
            names(B)[j] <- "phi";
            Bl[j] <- 0;
            Bu[j] <- 1;
            j <- j+1;
        }

        # Initials
        if(initialType=="optimal"){
            i <- 1;
            B[j] <- matVt[i,lagsModelMax];
            names(B)[j] <- "level";
            if(Etype=="A"){
                Bl[j] <- -Inf;
                Bu[j] <- Inf;
            }
            else{
                Bl[j] <- 0;
                Bu[j] <- Inf;
            }
            j <- j+1;
            i <- i+1;
            if(Ttype!="N"){
                B[j] <- matVt[i,lagsModelMax];
                names(B)[j] <- "trend";
                if(Ttype=="A"){
                    Bl[j] <- -Inf;
                    Bu[j] <- Inf;
                }
                else{
                    Bl[j] <- 0;
                    Bu[j] <- Inf;
                }
                j <- j+1;
                i <- i+1;
            }
            if(Stype!="N"){
                for(k in i:componentsNumber){
                    B[j+0:(lagsModel[k]-1)] <- matVt[k,(lagsModelMax-lagsModel[k])+1:lagsModel[k]];
                    names(B)[j+0:(lagsModel[k]-1)] <- paste0("seasonal",k,"_",1:lagsModel[k]);
                    if(Stype=="A"){
                        Bl[j+0:(lagsModel[k]-1)] <- -Inf;
                        Bu[j+0:(lagsModel[k]-1)] <- Inf;
                    }
                    else{
                        Bl[j+0:(lagsModel[k]-1)] <- 0;
                        Bu[j+0:(lagsModel[k]-1)] <- Inf;
                    }
                    j <- j+lagsModel[k];
                }
            }
        }

        # Initials of the xreg
        if(xregInitialsEstimate){
            B[j-1+1:xregNumber] <- matVt[componentsNumber+1:xregNumber,lagsModelMax];
            names(B)[j-1+1:xregNumber] <- rownames(matVt)[componentsNumber+1:xregNumber];
            if(Etype=="A"){
                Bl[j-1+1:xregNumber] <- -Inf;
                Bu[j-1+1:xregNumber] <- Inf;
            }
            else{
                Bl[j-1+1:xregNumber] <- -Inf;
                Bu[j-1+1:xregNumber] <- Inf;
            }
            j <- j+xregNumber;
        }

        # Add lambda if it is needed
        if(lambdaEstimate){
            B[j] <- 0.5;
            names(B)[j] <- "lambda";
            Bl[j] <- 1e-10;
            Bu[j] <- Inf;
        }

        return(list(B=B,Bl=Bl,Bu=Bu));
    }

    ##### Function returns scale parameter for the provided parameters #####
    scaler <- function(distribution, Etype, errors, yFitted, obsInSample, lambda){
        scale <- switch(distribution,
                        "dnorm"=,
                        "dt"=sqrt(sum(errors^2)/obsInSample),
                        "dlogis"=sqrt(sum(errors^2)/obsInSample * 3 / pi^2),
                        "dlaplace"=sum(abs(errors))/obsInSample,
                        "ds"=sum(sqrt(abs(errors))) / (obsInSample*2),
                        "dalaplace"=sum(errors*(lambda-(errors<=0)*1))/obsInSample,
                        "dlnorm"=switch(Etype,
                                        "A"=sqrt(sum(log(1+errors/yFitted)^2)/obsInSample),
                                        "M"=sqrt(sum(log(1+errors)^2)/obsInSample)),
                        "dinvgauss"=switch(Etype,
                                           "A"=sum((errors/yFitted)^2/(1+errors/yFitted))/obsInSample,
                                           "M"=sum((errors)^2/(1+errors))/obsInSample),
                        );
        return(scale);
    }

    ##### Cost Function for ETS #####
    CF <- function(B,
                   Etype, Ttype, Stype, yInSample,
                   ot, otLogical, occurrenceModel, obsInSample,
                   componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                   matVt, matWt, matF, vecG, componentsNumberSeasonal,
                   persistenceEstimate, phiEstimate, initialType,
                   xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                   xregNumber,
                   bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate){

        # Fill in the matrices
        mesElements <- filler(B,
                              Ttype, Stype, componentsNumber, lagsModel, lagsModelMax,
                              matVt, matWt, matF, vecG,
                              persistenceEstimate, phiEstimate, initialType,
                              xregInitialsEstimate, xregPersistenceEstimate, xregNumber);

        # If we estimate lambda, take it from the B vector
        if(lambdaEstimate){
            lambda[] <- B[length(B)];
        }

        # Check the bounds
        if(bounds=="usual"){
            if(any(mesElements$vecG>1) || any(mesElements$vecG<0)){
                return(1E+300);
            }
            # This is the restriction on the damping parameter
            if(mesElements$matF[2,2]>1 || mesElements$matF[2,2]<0){
                return(1E+300);
            }
        }
        else if(bounds=="admissible"){
            # We check the condition only for the last row of matWt
            eigenValues <- eigen(mesElements$matF - mesElements$vecG %*% mesElements$matWt[obsInSample,],
                                 only.values=TRUE)$values;
            if(any(abs(eigenValues)>1+1E-50)){
                return(abs(eigenValues)*1E+100);
            }
        }

        # Produce fitted values and errors
        mesFitted <- mesFitterWrap(mesElements$matVt, mesElements$matWt, mesElements$matF, mesElements$vecG,
                                   lagsModelAll, Etype, Ttype, Stype, componentsNumber, componentsNumberSeasonal,
                                   yInSample, ot, initialType=="backcasting");

        if(!multisteps){
            if(loss=="likelihood"){
                # Scale for different functions
                scale <- scaler(distribution, Etype, mesFitted$errors[otLogical],
                                mesFitted$yFitted[otLogical], obsInSample, lambda);

                # Calculate the likelihood
                CFValue <- -sum(switch(distribution,
                                       "dnorm"=switch(Etype,
                                                      "A"=dnorm(x=yInSample[otLogical], mean=mesFitted$yFitted[otLogical],
                                                                sd=scale, log=TRUE),
                                                      "M"=dnorm(x=yInSample[otLogical], mean=mesFitted$yFitted[otLogical],
                                                                sd=scale*mesFitted$yFitted[otLogical], log=TRUE)),
                                       "dlogis"=switch(Etype,
                                                       "A"=dlogis(x=yInSample[otLogical], location=mesFitted$yFitted[otLogical],
                                                                  scale=scale, log=TRUE),
                                                       "M"=dlogis(x=yInSample[otLogical], location=mesFitted$yFitted[otLogical],
                                                                  scale=scale*mesFitted$yFitted[otLogical], log=TRUE)),
                                       "dlaplace"=switch(Etype,
                                                         "A"=dlaplace(q=yInSample[otLogical], mu=mesFitted$yFitted[otLogical],
                                                                      scale=scale, log=TRUE),
                                                         "M"=dlaplace(q=yInSample[otLogical], mu=mesFitted$yFitted[otLogical],
                                                                      scale=scale*mesFitted$yFitted[otLogical], log=TRUE)),
                                       "dt"=switch(Etype,
                                                   "A"=dt(mesFitted$errors[otLogical], df=abs(lambda), log=TRUE),
                                                   "M"=dt(mesFitted$errors[otLogical]*mesFitted$yFitted[otLogical],
                                                          df=abs(lambda), log=TRUE)),
                                       "ds"=switch(Etype,
                                                   "A"=ds(q=yInSample[otLogical],mu=mesFitted$yFitted[otLogical],
                                                          scale=scale, log=TRUE),
                                                   "M"=ds(q=yInSample[otLogical],mu=mesFitted$yFitted[otLogical],
                                                          scale=scale*sqrt(mesFitted$yFitted[otLogical]), log=TRUE)),
                                       "dalaplace"=switch(Etype,
                                                          "A"=dalaplace(q=yInSample[otLogical], mu=mesFitted$yFitted[otLogical],
                                                                        scale=scale, alpha=lambda, log=TRUE),
                                                          "M"=dalaplace(q=yInSample[otLogical], mu=mesFitted$yFitted[otLogical],
                                                                        scale=scale*mesFitted$yFitted[otLogical], alpha=lambda, log=TRUE)),
                                       "dlnorm"=dlnorm(x=yInSample[otLogical], meanlog=log(mesFitted$yFitted[otLogical]),
                                                       sdlog=scale, log=TRUE),
                                       "dinvgauss"=dinvgauss(x=yInSample[otLogical], mean=mesFitted$yFitted[otLogical],
                                                             dispersion=scale/mesFitted$yFitted[otLogical], log=TRUE)));

                # Differential entropy for the logLik of occurrence model
                if(occurrenceModel){
                    CFValue <- CFValue + switch(distribution,
                                                "dnorm" =,
                                                "dlnorm" = obsZero*(log(sqrt(2*pi)*scale)+0.5),
                                                "dlogis" = obsZero*2,
                                                "dlaplace" =,
                                                "dalaplace" = obsZero*(1 + log(2*scale)),
                                                "dt" = obsZero*((scale+1)/2 *
                                                                    (digamma((scale+1)/2)-digamma(scale/2)) +
                                                                    log(sqrt(scale) * beta(scale/2,0.5))),
                                                "ds" = obsZero*(2 + 2*log(2*scale)),
                                                "dinvgauss" = obsZero*(0.5*(log(pi/2)+1+suppressWarnings(log(scale)))));
                }
            }
            else if(loss=="MSE"){
                CFValue <- mean(mesFitted$errors^2);
            }
            else if(loss=="MAE"){
                CFValue <- mean(abs(mesFitted$errors));
            }
            else if(loss=="HAM"){
                CFValue <- mean(sqrt(abs(mesFitted$errors)));
            }
            else if(any(loss==c("LASSO","RIDGE"))){
                ### All of this is needed in order to normalise level, trend, seasonal and xreg parameters
                # Define, how many elements to skip (we don't normalise smoothing parameters)
                if(xregPersistenceEstimate){
                    persistenceToSkip <- componentsNumber+xregNumber;
                }
                else{
                    persistenceToSkip <- componentsNumber;
                }
                j <- 1;
                if(phiEstimate){
                    j[] <- 2;
                }
                if(initialType=="optimal"){
                    # Standardise the level
                    if(Etype=="M"){
                        B[persistenceToSkip+j] <- log(B[persistenceToSkip+j] / mean(yInSample[1:lagsModelMax]));
                    }
                    else{
                        B[persistenceToSkip+j] <- B[persistenceToSkip+j] / mean(yInSample[1:lagsModelMax]);
                    }
                    # Change B values for the trend, so that it shrinks properly
                    if(Ttype=="M"){
                        j[] <- j+1;
                        B[persistenceToSkip+j] <- log(B[persistenceToSkip+j]);
                    }
                    else if(Ttype=="A"){
                        j[] <- j+1;
                        B[persistenceToSkip+j] <- B[persistenceToSkip+j]/mean(yInSample);
                    }
                    # Change B values for seasonality, so that it shrinks properly
                    if(Stype=="M"){
                        B[persistenceToSkip+j+1:(sum(lagsModel)-j)] <- log(B[persistenceToSkip+j+1:(sum(lagsModel)-j)]);
                    }
                    else if(Stype=="A"){
                        B[persistenceToSkip+j+1:(sum(lagsModel)-j)] <- B[persistenceToSkip+j+1:(sum(lagsModel)-j)]/mean(yInSample);
                    }

                    # Normalise parameters of xreg if they are additive. Otherwise leave - they will be small and close to zero
                    if(xregNumber>0 && Etype=="A"){
                        denominator <- tail(colMeans(matWt),xregNumber);
                        # If it is lower than 1, then we are probably dealing with (0, 1). No need to normalise
                        denominator[abs(denominator)<1] <- 1;
                        B[persistenceToSkip+sum(lagsModel)+c(1:xregNumber)] <- tail(B,xregNumber) / denominator;
                    }
                }

                CFValue <- switch(loss,
                                  "LASSO" = switch(Etype,
                                                   "A"=(1-lambda)* sqrt(sum(mesFitted$errors^2))/obsInSample +
                                                       lambda * sum(abs(B)),
                                                   "M"=(1-lambda)* sqrt(sum(log(1+mesFitted$errors)^2))/obsInSample +
                                                       lambda * sum(abs(B))),
                                  "RIDGE" = switch(Etype,
                                                   "A"=(1-lambda)* sqrt(sum(mesFitted$errors^2)) + lambda * sqrt(sum((B)^2)),
                                                   "M"=(1-lambda)* sqrt(sum(log(1+mesFitted$errors)^2)) + lambda * sqrt(sum((B)^2))));
            }
        }
        else{
            # Call for the Rcpp function to produce a matrix of multistep errors
            # The function should produce vector forecasts for the whole sample (Nikos' style)
            # loss==c("MSEh","TMSE","GTMSE","MSCE","MAEh","TMAE","GTMAE","MACE",
            #         "HAMh","THAM","GTHAM","CHAM","GPL",
            #         "aMSEh","aTMSE","aGTMSE","aMSCE","aGPL")
        }

        return(CFValue);
    }

    #### The function returns log-likelihood of the model ####
    logLikMES <- function(B,
                          Etype, Ttype, Stype, yInSample,
                          ot, otLogical, occurrenceModel, pFitted, obsInSample,
                          componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                          matVt, matWt, matF, vecG, componentsNumberSeasonal,
                          persistenceEstimate, phiEstimate, initialType,
                          xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                          xregNumber,
                          bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate){
        if(!multisteps){
            if(any(loss==c("LASSO","RIDGE"))){
                return(0);
            }
            else{
                distributionNew <- switch(loss,
                                          "MSE"=switch(Etype,"A"="dnorm","M"="dlnorm"),
                                          "MAE"="dlaplace",
                                          "HAM"="ds",
                                          distribution);
                logLikReturn <- -CF(B,
                                    Etype, Ttype, Stype, yInSample,
                                    ot, otLogical, occurrenceModel, obsInSample,
                                    componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                                    matVt, matWt, matF, vecG, componentsNumberSeasonal,
                                    persistenceEstimate, phiEstimate, initialType,
                                    xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                                    xregNumber,
                                    bounds, "likelihood", distributionNew, horizon, multisteps, lambda, lambdaEstimate);

                # If this is an occurrence model, add the probabilities
                if(occurrenceModel){
                    if(is.infinite(logLikReturn)){
                        logLikReturn[] <- 0;
                    }
                    if(any(c(1-pFitted[!otLogical]==0,pFitted[otLogical]==0))){
                        # return(-Inf);
                        ptNew <- pFitted[(pFitted!=0) & (pFitted!=1)];
                        otNew <- ot[(pFitted!=0) & (pFitted!=1)];

                        # Just return the original likelihood if the probability is weird
                        if(length(ptNew)==0){
                            return(logLikReturn);
                        }
                        else{
                            return(logLikReturn + sum(log(ptNew[otNew==1])) + sum(log(1-ptNew[otNew==0])));
                        }
                    }
                    else{
                        return(logLikReturn + sum(log(pFitted[otLogical])) + sum(log(1-pFitted[!otLogical])));
                    }
                }
                else{
                    return(logLikReturn);
                }
            }
        }
        else{
            # Use the predictive likelihoods from the GPL paper:
            # - Normal for MSEh, MSCE, GPL and their analytical counterparts
            # - Laplace for MAEh and MACE,
            # - S for HAMh and CHAM
        }
    }

    #### The function estimates the ETS model and returns B, logLik, nParam and CF(B) ####
    estimator <- function(Etype, Ttype, Stype, lags,
                          obsStates, obsInSample,
                          yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                          initialType, initialValue,
                          xregProvided, xregInitialsProvided, xregInitialsEstimate,
                          xregPersistence, xregPersistenceEstimate,
                          xregModel, xregData, xregNumber, xregNames,
                          ot, otLogical, occurrenceModel, pFitted,
                          bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate){

        # Create the basic variables
        mesArchitect <- architector(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType);
        list2env(mesArchitect, environment());

        # Create the matrices for the specific ETS model
        mesCreated <- creator(Etype, Ttype, Stype,
                              lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                              obsStates, obsInSample, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregProvided, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);

        if(is.null(B) && is.null(lb) && is.null(ub)){
            BValues <- initialiser(Etype, Ttype, Stype, componentsNumberSeasonal,
                                   componentsNumber, lagsModel, lagsModelMax, mesCreated$matVt,
                                   persistenceEstimate, phiEstimate, initialType,
                                   xregInitialsEstimate, xregPersistenceEstimate, xregNumber, lambdaEstimate);
            # Create the vector of initials for the optimisation
            B <- BValues$B;
            # lb <- BValues$Bl;
            # ub <- BValues$Bu;
        }

        # If the distribution is default, change it according to the error term
        if(loss=="likelihood" && distribution=="default"){
            distributionNew <- switch(Etype,
                                      "A"="dnorm",
                                      "M"="dinvgauss");
        }
        else{
            distributionNew <- distribution;
        }

        # Parameters are chosen to speed up the optimisation process and have decent accuracy
        res <- suppressWarnings(nloptr(B, CF, lb=lb, ub=ub,
                                       opts=list(algorithm=algorithm, xtol_rel=xtol_rel, maxeval=maxeval,
                                                 maxtime=maxtime, print_level=print_level),
                                       Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                                       ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, obsInSample=obsInSample,
                                       componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                                       matVt=mesCreated$matVt, matWt=mesCreated$matWt, matF=mesCreated$matF, vecG=mesCreated$vecG,
                                       componentsNumberSeasonal=componentsNumberSeasonal,
                                       persistenceEstimate=persistenceEstimate, phiEstimate=phiEstimate, initialType=initialType,
                                       xregProvided=xregProvided, xregInitialsEstimate=xregInitialsEstimate,
                                       xregPersistenceEstimate=xregPersistenceEstimate, xregNumber=xregNumber,
                                       bounds=bounds, loss=loss, distribution=distributionNew, horizon=horizon, multisteps=multisteps,
                                       lambda=lambda, lambdaEstimate=lambdaEstimate));

        ##### !!! Check the obtained parameters and the loss value and remove redundant parameters !!! #####
        # Cases to consider:
        # 1. Some smoothing parameters are zero or one;
        # 2. The cost function value is -Inf (due to no variability in the sample);

        # Prepare the values to return
        B[] <- res$solution;
        CFValue <- res$objective;
        # In case of likelihood, we typically have one more parameter to estimate - scale
        nParamEstimated <- length(B) + (loss=="likelihood");
        # Return a proper logLik class
        logLikMESValue <- structure(logLikMES(B,
                                              Etype, Ttype, Stype, yInSample,
                                              ot, otLogical, occurrenceModel, pFitted, obsInSample,
                                              componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                                              mesCreated$matVt, mesCreated$matWt, mesCreated$matF, mesCreated$vecG, componentsNumberSeasonal,
                                              persistenceEstimate, phiEstimate, initialType,
                                              xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                                              xregNumber,
                                              bounds, loss, distributionNew, horizon, multisteps, lambda, lambdaEstimate)
                                    ,nobs=obsInSample,df=nParamEstimated,class="logLik");

        return(list(B=B, CFValue=CFValue, nParamEstimated=nParamEstimated, logLikMESValue=logLikMESValue));
    }


    #### The function creates a pool of models and selects the best of them ####
    selector <- function(model, modelsPool, allowMultiplicative,
                         Etype, Ttype, Stype, damped, lags,
                         obsStates, obsInSample,
                         yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                         initialType, initialValue,
                         xregProvided, xregInitialsProvided, xregInitialsEstimate,
                         xregPersistence, xregPersistenceEstimate,
                         xregModel, xregData, xregNumber, xregNames,
                         ot, otLogical, occurrenceModel, pFitted, ICFunction,
                         bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate){

        # Check if the pool was provided. In case of "no", form the big and the small ones
        if(is.null(modelsPool)){
            # The variable saying that the pool was not provided.
            if(!silent){
                cat("Forming the pool of models based on... ");
            }

            # Define the whole pool of errors
            if(!allowMultiplicative){
                poolErrors <- c("A");
                poolTrends <- c("N","A","Ad");
                poolSeasonals <- c("N","A");
            }
            else{
                poolErrors <- c("A","M");
                poolTrends <- c("N","A","Ad","M","Md");
                poolSeasonals <- c("N","A","M");
            }

            # Some preparation variables
            # If Etype is not Z, then check on additive errors
            if(Etype!="Z"){
                poolErrors <- poolErrorsSmall <- Etype;
            }
            else{
                poolErrorsSmall <- "A";
            }

            # If Ttype is not Z, then create a pool with specified type
            if(Ttype!="Z"){
                if(Ttype=="X"){
                    poolTrendsSmall <- c("N","A");
                    poolTrends <- c("N","A","Ad");
                    checkTrend <- TRUE;
                }
                else if(Ttype=="Y"){
                    poolTrendsSmall <- c("N","M");
                    poolTrends <- c("N","M","Md");
                    checkTrend <- TRUE;
                }
                else{
                    if(damped){
                        poolTrends <- poolTrendsSmall <- paste0(Ttype,"d");
                    }
                    else{
                        poolTrends <- poolTrendsSmall <- Ttype;
                    }
                    checkTrend <- FALSE;
                }
            }
            else{
                poolTrendsSmall <- c("N","A");
                checkTrend <- TRUE;
            }

            # If Stype is not Z, then crete specific pools
            if(Stype!="Z"){
                if(Stype=="X"){
                    poolSeasonals <- poolSeasonalsSmall <- c("N","A");
                    checkSeasonal <- TRUE;
                }
                else if(Stype=="Y"){
                    poolSeasonalsSmall <- c("N","M");
                    poolSeasonals <- c("N","M");
                    checkSeasonal <- TRUE;
                }
                else{
                    poolSeasonalsSmall <- Stype;
                    poolSeasonals <- Stype;
                    checkSeasonal <- FALSE;
                }
            }
            else{
                poolSeasonalsSmall <- c("N","A","M");
                checkSeasonal <- TRUE;
            }

            # If ZZZ, then the vector is: "ANN" "ANA" "ANM" "AAN" "AAA" "AAM"
            # Otherwise id depends on the provided restrictions
            poolSmall <- paste0(rep(poolErrorsSmall,length(poolTrendsSmall)*length(poolSeasonalsSmall)),
                                 rep(poolTrendsSmall,each=length(poolSeasonalsSmall)),
                                 rep(poolSeasonalsSmall,length(poolTrendsSmall)));
            modelsTested <- NULL;
            modelCurrent <- NA;

            # Counter + checks for the components
            j <- 1;
            i <- 0;
            check <- TRUE;
            besti <- bestj <- 1;
            results <- vector("list",length(poolSmall));

            #### Branch and bound is here ####
            while(check){
                i <- i + 1;
                modelCurrent[] <- poolSmall[j];
                if(!silent){
                    cat(paste0(modelCurrent,", "));
                }
                Etype[] <- substring(modelCurrent,1,1);
                Ttype[] <- substring(modelCurrent,2,2);
                if(nchar(modelCurrent)==4){
                    phi[] <- 0.95;
                    phiEstimate[] <- TRUE;
                    Stype[] <- substring(modelCurrent,4,4);
                }
                else{
                    phi <- 1;
                    phiEstimate[] <- FALSE;
                    Stype[] <- substring(modelCurrent,3,3);
                }

                results[[i]] <- estimator(Etype, Ttype, Stype, lags,
                                          obsStates, obsInSample,
                                          yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                          initialType, initialValue,
                                          xregProvided, xregInitialsProvided, xregInitialsEstimate,
                                          xregPersistence, xregPersistenceEstimate,
                                          xregModel, xregData, xregNumber, xregNames,
                                          ot, otLogical, occurrenceModel, pFitted,
                                          bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);
                results[[i]]$IC <- ICFunction(results[[i]]$logLikMESValue);
                results[[i]]$Etype <- Etype;
                results[[i]]$Ttype <- Ttype;
                results[[i]]$Stype <- Stype;
                results[[i]]$phiEstimate <- phiEstimate;
                results[[i]]$model <- modelCurrent;

                modelsTested <- c(modelsTested,modelCurrent);

                if(j>1){
                    # If the first is better than the second, then choose first
                    if(results[[besti]]$IC <= results[[i]]$IC){
                        # If Ttype is the same, then we check seasonality
                        if(substring(modelCurrent,2,2)==substring(poolSmall[bestj],2,2)){
                            poolSeasonals <- results[[besti]]$Stype;
                            checkSeasonal <- FALSE;
                            # j[] <- j+1;
                            j <- which(poolSmall!=poolSmall[bestj] &
                                           substring(poolSmall,nchar(poolSmall),nchar(poolSmall))==poolSeasonals);
                        }
                        # Otherwise we checked trend
                        else{
                            poolTrends <- results[[bestj]]$Ttype;
                            checkTrend[] <- FALSE;
                        }
                    }
                    else{
                        if(substring(modelCurrent,2,2) == substring(poolSmall[besti],2,2)){
                            poolSeasonals <- poolSeasonals[poolSeasonals!=results[[besti]]$Stype];
                            if(length(poolSeasonals)>1){
                                # Select another seasonal model, that is not from the previous iteration and not the current one
                                bestj[] <- j;
                                besti[] <- i;
                                # j[] <- 3;
                                j <- 3;
                            }
                            else{
                                bestj[] <- j;
                                besti[] <- i;
                                j <- which(substring(poolSmall,nchar(poolSmall),nchar(poolSmall))==poolSeasonals &
                                               substring(poolSmall,2,2)!=substring(modelCurrent,2,2));
                                checkSeasonal[] <- FALSE;
                            }
                        }
                        else{
                            poolTrends <- poolTrends[poolTrends!=results[[bestj]]$Ttype];
                            besti[] <- i;
                            bestj[] <- j;
                            checkTrend[] <- FALSE;
                        }
                    }

                    if(all(!c(checkTrend,checkSeasonal))){
                        check[] <- FALSE;
                    }
                }
                else{
                    j <- 2;
                }

                if(j>=length(poolSmall)){
                    check[] <- FALSE;
                }
            }

            # Prepare a bigger pool based on the small one
            modelsPool <- unique(c(modelsTested,
                                   paste0(rep(poolErrors,each=length(poolTrends)*length(poolSeasonals)),
                                          poolTrends,
                                          rep(poolSeasonals,each=length(poolTrends)))));
            j <- length(modelsTested);
        }
        else{
            j <- 0;
            results <- vector("list",length(modelsPool));
        }
        modelsNumber <- length(modelsPool);

        #### Run the full pool of models ####
        if(!silent){
            cat("Estimation progress:    ");
        }
        # Start loop of models
        while(j < modelsNumber){
            j <- j + 1;
            if(!silent){
                if(j==1){
                    cat("\b");
                }
                cat(paste0(rep("\b",nchar(round((j-1)/modelsNumber,2)*100)+1),collapse=""));
                cat(paste0(round(j/modelsNumber,2)*100,"%"));
            }

            modelCurrent <- modelsPool[j];
            Etype <- substring(modelCurrent,1,1);
            Ttype <- substring(modelCurrent,2,2);
            if(nchar(modelCurrent)==4){
                phi <- 0.95;
                Stype <- substring(modelCurrent,4,4);
            }
            else{
                phi <- 1;
                Stype <- substring(modelCurrent,3,3);
            }

            results[[j]] <- estimator(Etype, Ttype, Stype, lags,
                                      obsStates, obsInSample,
                                      yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                      initialType, initialValue,
                                      xregProvided, xregInitialsProvided, xregInitialsEstimate,
                                      xregPersistence, xregPersistenceEstimate,
                                      xregModel, xregData, xregNumber, xregNames,
                                      ot, otLogical, occurrenceModel, pFitted,
                                      bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);
            results[[j]]$IC <- ICFunction(results[[j]]$logLikMESValue);
            results[[j]]$Etype <- Etype;
            results[[j]]$Ttype <- Ttype;
            results[[j]]$Stype <- Stype;
            results[[j]]$phiEstimate <- phiEstimate;
            results[[j]]$model <- modelCurrent;
        }

        if(!silent){
            cat("... Done! \n");
        }

        # Extract ICs and find the best
        icSelection <- vector("numeric",modelsNumber);
        for(i in 1:modelsNumber){
            icSelection[i] <- results[[i]]$IC;
        }
        names(icSelection) <- modelsPool;

        icSelection[is.nan(icSelection)] <- 1E100;

        return(list(results=results,icSelection=icSelection));
    }

    ##### !!!! This function will use residuals in order to determine the needed xreg !!!! #####
    # XregSelector <- function(listToReturn){
    # }

    ##### Function prepares all the matrices and vectors for return #####
    preparator <- function(B, Etype, Ttype, Stype,
                           lagsModel, lagsModelMax, lagsModelAll,
                           componentsNumber, componentsNumberSeasonal,
                           xregNumber, distribution, loss,
                           persistenceEstimate, phiEstimate, lambdaEstimate, initialType,
                           xregInitialsEstimate, xregPersistenceEstimate,
                           matVt, matWt, matF, vecG,
                           occurrenceModel, ot, oesModel,
                           parametersNumber, CFValue){

        # Fill in the matrices
        mesElements <- filler(B,
                              Ttype, Stype, componentsNumber, lagsModel, lagsModelMax,
                              matVt, matWt, matF, vecG,
                              persistenceEstimate, phiEstimate, initialType,
                              xregInitialsEstimate, xregPersistenceEstimate, xregNumber);
        list2env(mesElements, environment());

        # Write down lambda
        if(lambdaEstimate){
            lambda[] <- tail(B,1);
        }

        # Write down phi
        if(phiEstimate){
            phi[] <- B[names(B)=="phi"];
        }

        # Fit the model to the data
        mesFitted <- mesFitterWrap(matVt, matWt, matF, vecG,
                                   lagsModelAll, Etype, Ttype, Stype, componentsNumber, componentsNumberSeasonal,
                                   yInSample, ot, initialType=="backcasting");

        errors <- mesFitted$errors;
        yFitted <- mesFitted$yFitted;
        if(occurrenceModel){
            yFitted[] <- yFitted * pFitted;
        }

        matVt[] <- mesFitted$matVt;
        if(initialType=="backcasting"){
            matVt <- matVt[,1:(obsInSample+lagsModelMax), drop=FALSE];
        }

        if(horizon>0){
            yForecast <- ts(rep(NA, horizon), start=time(y)[obsInSample]+deltat(y), frequency=frequency(y));
            mesForecast <- mesForecasterWrap(matVt[,obsInSample+(1:lagsModelMax),drop=FALSE], tail(matWt,horizon), matF, vecG,
                                             lagsModelAll, Etype, Ttype, Stype,
                                             componentsNumber, componentsNumberSeasonal, horizon);
            yForecast[] <- mesForecast$yForecast;
            if(occurrenceModel){
                yForecast[] <- yForecast * forecast(oesModel, h=h)$mean;
            }
        }
        else{
            yForecast <- ts(NA, start=time(y)[obsInSample]+deltat(y), frequency=frequency(y));
        }

        # If the distribution is default, change it according to the error term
        if(loss=="likelihood" && distribution=="default"){
            distribution[] <- switch(Etype,
                                     "A"="dnorm",
                                     "M"="dinvgauss");
        }

        if(initialType=="optimal"){
            initialValue <- vector("numeric", sum(lagsModelAll));
            j <- 0;
            for(i in 1:length(lagsModelAll)){
                initialValue[j+1:lagsModelAll[i]] <- matVt[i,1:lagsModelAll[i]];
                j <- j + lagsModelAll[i];
            }
        }

        if(persistenceEstimate){
            persistence <- vecG[1:componentsNumber,];
            names(persistence) <- rownames(vecG)[1:componentsNumber];
        }

        if(xregPersistenceEstimate){
            xregPersistence <- vecG[-c(1:componentsNumber),];
            names(xregPersistence) <- rownames(vecG)[-c(1:componentsNumber)];
        }

        scale <- scaler(distribution, Etype, errors[otLogical], yFitted[otLogical], obsInSample, lambda);
        yFitted <- ts(yFitted,start=start(y), frequency=frequency(y));

        return(list(model=NA, timeElapsed=NA,
                    y=NA, holdout=NA, fitted=yFitted, residuals=errors,
                    forecast=yForecast, states=ts(t(matVt), start=(time(y)[1] - deltat(y)*lagsModelMax),
                                                  frequency=frequency(y)),
                    persistence=persistence, phi=phi, transition=matF,
                    measurement=matWt, initialType=initialType, initial=initialValue,
                    nParam=parametersNumber, occurrence=oesModel, xreg=xreg,
                    xregInitial=xregInitial, xregPersistence=xregPersistence,
                    loss=loss, lossValue=CFValue, logLik=logLikMESValue, distribution=distribution,
                    scale=scale, lambda=lambda, B=B, lags=lagsModel, FI=FI));
    }

    #### Deal with occurrence model ####
    if(occurrenceModel && !occurrenceModelProvided){
        oesModel <- suppressWarnings(oes(yInSample, model=model, occurrence=occurrence, ic=ic, h=horizon,
                                         holdout=FALSE, bounds="usual", xreg=xreg, xregDo=xregDo));
        pFitted[] <- fitted(oesModel);
        parametersNumber[1,3] <- nparam(oesModel);
        # This should not happen, but just in case...
        if(oesModel$occurrence=="n"){
            occurrence <- "n";
            otLogical <- rep(TRUE,obsInSample);
            occurrenceModel <- FALSE;
            ot <- matrix(otLogical*1,ncol=1);
            obsNonzero <- sum(ot);
            obsZero <- obsInSample - obsNonzero;
            Etype[] <- switch(Etype,
                              "M"="A",
                              "Y"=,
                              "Z"="X",
                              Etype);
            Ttype[] <- switch(Ttype,
                              "M"="A",
                              "Y"=,
                              "Z"="X",
                              Ttype);
            Stype[] <- switch(Stype,
                              "M"="A",
                              "Y"=,
                              "Z"="X",
                              Stype);
        }
    }
    else if(occurrenceModel && occurrenceModelProvided){
        parametersNumber[2,3] <- nparam(oesModel);
    }

    ##### Estimate the specified model #####
    if(modelDo=="estimate"){
        # Estimate the parameters of the demand sizes model
        mesEstimated <- estimator(Etype, Ttype, Stype, lags,
                                 obsStates, obsInSample,
                                 yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                 initialType, initialValue,
                                 xregProvided, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames,
                                 ot, otLogical, occurrenceModel, pFitted,
                                 bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);
        list2env(mesEstimated, environment());

        #### This part is needed in order for the filler to do its job later on
        # Create the basic variables based on the estimated model
        mesArchitect <- architector(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType);
        list2env(mesArchitect, environment());

        # Create the matrices for the specific ETS model
        mesCreated <- creator(Etype, Ttype, Stype,
                              lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                              obsStates, obsInSample, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregProvided, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);
        list2env(mesCreated, environment());

        ####!!! If the occurrence is auto, then compare this with the model with no occurrence !!!####

        parametersNumber[1,1] <- (sum(lagsModel)*(initialType=="optimal") + phiEstimate +
                                      componentsNumber*persistenceEstimate + xregNumber*xregInitialsEstimate +
                                      xregNumber*xregPersistenceEstimate + 1);
        if(xregProvided){
            parametersNumber[1,2] <- xregNumber*xregInitialsEstimate + xregNumber*xregPersistenceEstimate;
        }
        parametersNumber[1,4] <- sum(parametersNumber[1,1:3]);
        parametersNumber[2,4] <- sum(parametersNumber[2,1:3]);
    }
    #### Selection of the best model ####
    else if(modelDo=="select"){
        mesSelected <-  selector(model, modelsPool, allowMultiplicative,
                                 Etype, Ttype, Stype, damped, lags,
                                 obsStates, obsInSample,
                                 yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                 initialType, initialValue,
                                 xregProvided, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames,
                                 ot, otLogical, occurrenceModel, pFitted, ICFunction,
                                 bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);

        icSelection <- mesSelected$icSelection;
        # Take the parameters of the best model
        list2env(mesSelected$results[[which.min(icSelection)[1]]], environment());

        #### This part is needed in order for the filler to do its job later on
        # Create the basic variables based on the estimated model
        mesArchitect <- architector(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType);
        list2env(mesArchitect, environment());

        # Create the matrices for the specific ETS model
        mesCreated <- creator(Etype, Ttype, Stype,
                              lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                              obsStates, obsInSample, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregProvided, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);
        list2env(mesCreated, environment());

        parametersNumber[1,1] <- (sum(lagsModel)*(initialType=="optimal") + phiEstimate +
                                      componentsNumber*persistenceEstimate + xregNumber*xregInitialsEstimate +
                                      xregNumber*xregPersistenceEstimate + 1);
        if(xregProvided){
            parametersNumber[1,2] <- xregNumber*xregInitialsEstimate + xregNumber*xregPersistenceEstimate;
        }
        parametersNumber[1,4] <- sum(parametersNumber[1,1:3]);
        parametersNumber[2,4] <- sum(parametersNumber[2,1:3]);
    }
    #### Combination of models ####
    else if(modelDo=="combine"){
        modelOriginal <- model;
        # If the pool is not provided, then create one
        if(is.null(modelsPool)){
            # Define the whole pool of errors
            if(!allowMultiplicative){
                poolErrors <- c("A");
                poolTrends <- c("N","A","Ad");
                poolSeasonals <- c("N","A");
            }
            else{
                poolErrors <- c("A","M");
                poolTrends <- c("N","A","Ad","M","Md");
                poolSeasonals <- c("N","A","M");
            }

            # Some preparation variables
            # If Etype is not Z, then check on additive errors
            if(Etype!="Z"){
                poolErrors <- switch(Etype,
                                     "N"="N",
                                     "A"=,
                                     "X"="A",
                                     "M"=,
                                     "Y"="M");
            }

            # If Ttype is not Z, then create a pool with specified type
            if(Ttype!="Z"){
                poolTrends <- switch(Ttype,
                                     "N"="N",
                                     "A"=ifelse(damped,"Ad","A"),
                                     "M"=ifelse(damped,"Md","M"),
                                     "X"=c("N","A","Ad"),
                                     "Y"=c("N","M","Md"));
            }

            # If Stype is not Z, then crete specific pools
            if(Stype!="Z"){
                poolSeasonals <- switch(Stype,
                                     "N"="N",
                                     "A"="A",
                                     "X"=c("N","A"),
                                     "M"="M",
                                     "Y"=c("N","M"));
            }

            modelsPool <- paste0(rep(poolErrors,length(poolTrends)*length(poolSeasonals)),
                                 rep(poolTrends,each=length(poolSeasonals)),
                                 rep(poolSeasonals,length(poolTrends)));
        }

        mesSelected <-  selector(model, modelsPool, allowMultiplicative,
                                 Etype, Ttype, Stype, damped, lags,
                                 obsStates, obsInSample,
                                 yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                 initialType, initialValue,
                                 xregProvided, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames,
                                 ot, otLogical, occurrenceModel, pFitted, ICFunction,
                                 bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);

        icSelection <- mesSelected$icSelection;

        icBest <- min(icSelection);
        mesSelected$icWeights  <- (exp(-0.5*(icSelection-icBest)) /
                                       sum(exp(-0.5*(icSelection-icBest))));

        for(i in 1:length(mesSelected$results)){
            # Take the parameters of the best model
            list2env(mesSelected$results[[i]], environment());

            #### This part is needed in order for the filler to do its job later on
            # Create the basic variables based on the estimated model
            mesArchitect <- architector(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType);
            list2env(mesArchitect, environment());

            mesSelected$results[[i]]$lagsModel <- mesArchitect$lagsModel;
            mesSelected$results[[i]]$lagsModelAll <- mesArchitect$lagsModelAll;
            mesSelected$results[[i]]$lagsModelMax <- mesArchitect$lagsModelMax;
            mesSelected$results[[i]]$lagsLength <- mesArchitect$lagsLength;
            mesSelected$results[[i]]$componentsNumber <- mesArchitect$componentsNumber;
            mesSelected$results[[i]]$componentsNumberSeasonal <- mesArchitect$componentsNumberSeasonal;
            mesSelected$results[[i]]$componentsNames <- mesArchitect$componentsNames;

            # Create the matrices for the specific ETS model
            mesCreated <- creator(Etype, Ttype, Stype,
                                  lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                                  obsStates, obsInSample, componentsNumber, componentsNumberSeasonal,
                                  componentsNames, otLogical,
                                  yInSample, persistence, persistenceEstimate, phi,
                                  initialValue, initialType,
                                  xregProvided, xregInitialsProvided, xregPersistence,
                                  xregModel, xregData, xregNumber, xregNames);

            mesSelected$results[[i]]$matVt <- mesCreated$matVt;
            mesSelected$results[[i]]$matWt <- mesCreated$matWt;
            mesSelected$results[[i]]$matF <- mesCreated$matF;
            mesSelected$results[[i]]$vecG <- mesCreated$vecG;

            parametersNumber[1,1] <- (sum(lagsModel)*(initialType=="optimal") + phiEstimate +
                                          componentsNumber*persistenceEstimate + xregNumber*xregInitialsEstimate +
                                          xregNumber*xregPersistenceEstimate + 1);
            if(xregProvided){
                parametersNumber[1,2] <- xregNumber*xregInitialsEstimate + xregNumber*xregPersistenceEstimate;
            }
            parametersNumber[1,4] <- sum(parametersNumber[1,1:3]);
            parametersNumber[2,4] <- sum(parametersNumber[2,1:3]);

            mesSelected$results[[i]]$parametersNumber <- parametersNumber;
        }
    }
    #### Use the provided model ####
    else if(modelDo=="use"){
        # If the distribution is default, change it according to the error term
        if(loss=="likelihood" && distribution=="default"){
            distributionNew <- switch(Etype,
                                      "A"="dnorm",
                                      "M"="dinvgauss");
        }
        else{
            distributionNew <- distribution;
        }

        # Create the basic variables
        mesArchitect <- architector(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType);
        list2env(mesArchitect, environment());

        # Create the matrices for the specific ETS model
        mesCreated <- creator(Etype, Ttype, Stype,
                              lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                              obsStates, obsInSample, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregProvided, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);
        list2env(mesCreated, environment());

        CFValue <- CF(B=0, Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                      ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, obsInSample=obsInSample,
                      componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                      matVt=mesCreated$matVt, matWt=mesCreated$matWt, matF=mesCreated$matF, vecG=mesCreated$vecG,
                      componentsNumberSeasonal=componentsNumberSeasonal,
                      persistenceEstimate=persistenceEstimate, phiEstimate=phiEstimate, initialType=initialType,
                      xregProvided=xregProvided, xregInitialsEstimate=xregInitialsEstimate,
                      xregPersistenceEstimate=xregPersistenceEstimate, xregNumber=xregNumber,
                      bounds=bounds, loss=loss, distribution=distributionNew, horizon=horizon, multisteps=multisteps,
                      lambda=lambda, lambdaEstimate=lambdaEstimate);

        parametersNumber[1,1] <- parametersNumber[1,4] <- 1;
        logLikMESValue <- structure(logLikMES(B,
                                              Etype, Ttype, Stype, yInSample,
                                              ot, otLogical, occurrenceModel, pFitted, obsInSample,
                                              componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                                              mesCreated$matVt, mesCreated$matWt, mesCreated$matF, mesCreated$vecG, componentsNumberSeasonal,
                                              persistenceEstimate, phiEstimate, initialType,
                                              xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                                              xregNumber,
                                              bounds, loss, distributionNew, horizon, multisteps, lambda, lambdaEstimate)
                                    ,nobs=obsInSample,df=parametersNumber[1,4],class="logLik")

        # If Fisher Information is required, do that analytically
        if(FI){
            # Define parameters just for FI calculation
            if(initialType=="provided" && any(names(B)=="level")){
                initialTypeFI <- "optimal";
            }
            else{
                initialTypeFI <- initialType;
            }

            #### This will not work in cases, when initials for xreg were provided from the beginning!
            if(xregProvided && initialTypeFI=="optimal"){
                xregInitialsEstimateFI <- TRUE;
            }
            else{
                xregInitialsEstimateFI <- FALSE;
            }
            # If smoothing parmaeters were estimated, then alpha should be in the list
            if(any(names(B)=="alpha")){
                persistenceEstimateFI <- TRUE;
            }
            else{
                persistenceEstimateFI <- FALSE;
            }
            if(any(names(B)=="phi")){
                phiEstimateFI <- TRUE;
            }
            else{
                phiEstimateFI <- FALSE;
            }
            if(any(names(B)=="lambda")){
                lambdaEstimateFI <- TRUE;
            }
            else{
                lambdaEstimateFI <- FALSE;
            }
            if(any(substr(names(B),1,5)=="delta")){
                xregPersistenceEstimateFI <- TRUE;
            }
            else{
                xregPersistenceEstimateFI <- FALSE;
            }

            FI <- hessian(CF, B, Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                          ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, obsInSample=obsInSample,
                          componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                          matVt=matVt, matWt=matWt, matF=matF, vecG=vecG,
                          componentsNumberSeasonal=componentsNumberSeasonal,
                          persistenceEstimate=persistenceEstimateFI, phiEstimate=phiEstimateFI, initialType=initialTypeFI,
                          xregProvided=xregProvided, xregInitialsEstimate=xregInitialsEstimateFI,
                          xregPersistenceEstimate=xregPersistenceEstimateFI, xregNumber=xregNumber,
                          bounds=bounds, loss=loss, distribution=distribution, horizon=horizon, multisteps=multisteps,
                          lambda=lambda, lambdaEstimate=lambdaEstimateFI);
            colnames(FI) <- names(B);
            rownames(FI) <- names(B);
        }
        else{
            FI <- NULL;
        }
    }

    # Transform everything into ts
    yInSample <- ts(yInSample,start=start(y), frequency=frequency(y));
    if(holdout){
        yHoldout <- ts(yHoldout, start=time(y)[obsInSample]+deltat(y), frequency=frequency(y));
    }

    #### Prepare the return if we didn't combine anything ####
    if(modelDo!="combine"){
        modelReturned <- preparator(B, Etype, Ttype, Stype,
                                    lagsModel, lagsModelMax, lagsModelAll,
                                    componentsNumber, componentsNumberSeasonal,
                                    xregNumber, distribution, loss,
                                    persistenceEstimate, phiEstimate, lambdaEstimate, initialType,
                                    xregInitialsEstimate, xregPersistenceEstimate,
                                    matVt, matWt, matF, vecG,
                                    occurrenceModel, ot, oesModel,
                                    parametersNumber, CFValue);

        # Prepare the name of the model
        if(xregExist){
            modelName <- "ETSX";
        }
        else{
            modelName <- "ETS";
        }
        modelName <- paste0(modelName,"(",model,")");
        if(all(occurrence!=c("n","none"))){
            modelName <- paste0("i",modelName);
        }
        if(componentsNumberSeasonal>1){
            modelName <- paste0(modelName,"[",paste0(lags[lags!=1], collapse=", "),"]");
        }

        modelReturned$model <- modelName;
        modelReturned$timeElapsed <- Sys.time()-startTime;
        modelReturned$y <- yInSample;
        modelReturned$holdout <- yHoldout;

        class(modelReturned) <- c("mes","smooth");
    }
    else{
        modelReturned <- list(models=vector("list",length(mesSelected$results)));
        yFittedCombined <- rep(0,obsInSample);
        if(h>0){
            yForecastCombined <- rep(0,h);
        }
        else{
            yForecastCombined <- NA;
        }
        parametersNumberOverall <- parametersNumber;

        for(i in 1:length(mesSelected$results)){
            list2env(mesSelected$results[[i]], environment());
            modelReturned$models[[i]] <- preparator(B, Etype, Ttype, Stype,
                                                   lagsModel, lagsModelMax, lagsModelAll,
                                                   componentsNumber, componentsNumberSeasonal,
                                                   xregNumber, distribution, loss,
                                                   persistenceEstimate, phiEstimate, lambdaEstimate, initialType,
                                                   xregInitialsEstimate, xregPersistenceEstimate,
                                                   matVt, matWt, matF, vecG,
                                                   occurrenceModel, ot, oesModel,
                                                   parametersNumber, CFValue);
            yFittedCombined[] <- yFittedCombined + modelReturned$models[[i]]$fitted * mesSelected$icWeights[i];
            if(h>0){
                yForecastCombined[] <- yForecastCombined + modelReturned$models[[i]]$forecast * mesSelected$icWeights[i];
            }

            # Prepare the name of the model
            if(xregExist){
                modelName <- "ETSX";
            }
            else{
                modelName <- "ETS";
            }
            modelName <- paste0(modelName,"(",model,")");
            if(all(occurrence!=c("n","none"))){
                modelName <- paste0("i",modelName);
            }
            if(componentsNumberSeasonal>1){
                modelName <- paste0(modelName,"[",paste0(lags[lags!=1], collapse=", "),"]");
            }
            modelReturned$models[[i]]$model <- modelName;
            modelReturned$models[[i]]$timeElapsed <- Sys.time()-startTime;
            parametersNumberOverall[1,1] <- parametersNumber[1,1] * mesSelected$icWeights[i];
            modelReturned$models[[i]]$y <- yInSample;

            class(modelReturned$models[[i]]) <- c("mes","smooth");
        }

        # Record the original name of the model.
        model[] <- modelOriginal;
        # Prepare the name of the model
        if(xregExist){
            modelName <- "ETSX";
        }
        else{
            modelName <- "ETS";
        }
        modelName <- paste0(modelName,"(",model,")");
        if(all(occurrence!=c("n","none"))){
            modelName <- paste0("i",modelName);
        }
        if(componentsNumberSeasonal>1){
            modelName <- paste0(modelName,"[",paste0(lags[lags!=1], collapse=", "),"]");
        }
        modelReturned$model <- modelName;
        modelReturned$timeElapsed <- Sys.time()-startTime;
        modelReturned$holdout <- yHoldout;
        modelReturned$y <- yInSample;
        modelReturned$fitted <- ts(yFittedCombined,start=start(y), frequency=frequency(y));
        if(h>0){
            modelReturned$forecast <- ts(yForecastCombined,start=time(y)[obsInSample]+deltat(y), frequency=frequency(y));
        }
        else{
            modelReturned$forecast <- yForecastCombined;
        }
        parametersNumberOverall[1,4] <- sum(parametersNumberOverall[1,1:4]);
        modelReturned$nParam <- parametersNumberOverall;
        modelReturned$ICw <- mesSelected$icWeights;
        class(modelReturned) <- c("mesCombined","mes","smooth");
    }

    # model <- structure(list(model=modelName, timeElapsed=Sys.time()-startTime,
    #                         y=yInSample, holdout=yHoldout, fitted=yFitted, residuals=errors,
    #                         forecast=yForecast, states=ts(t(matVt), start=(time(y)[1] - deltat(y)*lagsModelMax),
    #                                                       frequency=frequency(y)),
    #                         persistence=persistence, phi=phi, transition=matF,
    #                         measurement=matWt, initialType=initialType, initial=initialValue,
    #                         nParam=parametersNumber, occurrence=oesModel, xreg=xreg,
    #                         xregInitial=xregInitial, xregPersistence=xregPersistence,
    #                         loss=loss, lossValue=CFValue, logLik=logLikMESValue, distribution=distribution,
    #                         scale=scale, lambda=lambda, B=B, lags=lagsModel, FI=FI),
    #                    class=c("mes","smooth"));

    if(!silent){
        plot(modelReturned, 1);
    }

    return(modelReturned);
}

#### Methods for mes ####
confint.mes <- function(object, parm, level=0.95, ...){
    mesVcov <- vcov(object);
    mesSD <- sqrt(abs(diag(mesVcov)));
    # mesCoef <- coef(object);
    mesCoefBounds <- matrix(0,length(mesSD),2);
    mesCoefBounds[,1] <- qnorm((1-level)/2, 0, mesSD);
    mesCoefBounds[,2] <- qnorm((1+level)/2, 0, mesSD);
    mesReturn <- cbind(mesSD,mesCoefBounds);
    colnames(mesReturn) <- c("S.E.",
                             paste0((1-level)/2*100,"%"), paste0((1+level)/2*100,"%"));
    # colnames(mesReturn) <- c("Estimate", "Std. Error",
    #                          paste0("Lower ",(1-level)/2*100,"%"), paste0("Upper ",(1+level)/2*100,"%"));

    return(mesReturn);
}

#' @export
coef.mes <- function(object, ...){
    return(object$B);
}
#
# covar.mes <- function(object, type=c("analytical","empirical","simulated"), ...){
#     h <- length(object$holdout);
#     lagsModel <- lags(object);
#     s2 <- sigma(object)^2;
#     persistence <- matrix(object$persistence,length(object$persistence),1);
#     transition <- object$transition;
#     measurement <- object$measurement;
#
#     covarMat <- covarAnal(lagsModel, h, measurement, transition, persistence, s2);
# }

#' @export
lags.mes <- function(object, ...){
    return(object$lags);
}

#' @export
residuals.mes <- function(object, ...){
    return(switch(object$distribution,
                  "dnorm"=,
                  "dlogis"=,
                  "dlaplace"=,
                  "dt"=,
                  "ds"=,
                  "dalaplace"=object$residuals,
                  "dlnorm"=,
                  "dinvgauss"=switch(errorType(object),
                                     "A"=1+object$residuals/fitted(object),
                                     "M"=1+object$residuals)));
}

#' @importFrom stats rstandard
#' @export
rstandard.mes <- function(model, ...){
    obs <- nobs(model);
    df <- obs - nparam(model);
    errors <- residuals(model);
    # If this is an occurrence model, then only modify the non-zero obs
    if(is.oes(model$occurrence)){
        residsToGo <- which(actuals(model$occurrence)!=0);
    }
    else{
        residsToGo <- c(1:obs);
    }
    if(any(model$distribution==c("dt","dnorm"))){
        return((errors - mean(errors[residsToGo])) / sqrt(model$scale^2 * obs / df));
    }
    else if(model$distribution=="ds"){
        return((errors - mean(errors[residsToGo])) / (model$scale * obs / df)^2);
    }
    else if(model$distribution=="dinvgauss"){
        return(errors / mean(errors[residsToGo]));
    }
    else if(model$distribution=="dlnorm"){
        errors[] <- log(errors);
        return(exp((errors - mean(errors[residsToGo])) / sqrt(model$scale^2 * obs / df)));
    }
    else{
        return(errors / model$scale * obs / df);
    }
}

#' @importFrom stats rstudent
#' @export
rstudent.mes <- function(model, ...){
    obs <- nobs(model);
    df <- obs - nparam(model) - 1;
    rstudentised <- errors <- residuals(model);
    # If this is an occurrence model, then only modify the non-zero obs
    if(is.alm(model$occurrence)){
        residsToGo <- which(actuals(model$occurrence)!=0);
    }
    else{
        residsToGo <- c(1:obs);
    }
    if(any(model$distribution==c("dt","dnorm"))){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / sqrt(sum(errors[-i]^2) / df);
        }
    }
    else if(model$distribution=="ds"){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sum(sqrt(abs(errors[-i]))) / (2*df))^2;
        }
    }
    else if(model$distribution=="dlaplace"){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sum(abs(errors[-i])) / df);
        }
    }
    else if(model$distribution=="dalaplace"){
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sum(errors[-i] * (model$lambda - (errors[-i]<=0)*1)) / df);
        }
    }
    else if(model$distribution=="dlogis"){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sqrt(sum(errors[-i]^2) / df) * sqrt(3) / pi);
        }
    }
    else if(model$distribution=="dlnorm"){
        errors[] <- log(errors) - mean(log(errors));
        for(i in residsToGo){
            rstudentised[i] <- exp(errors[i] / sqrt(sum(errors[-i]^2) / df));
        }
    }
    else if(model$distribution=="dinvgauss"){
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / mean(errors[residsToGo][-i]);
        }
    }
    else{
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / sqrt(sum(errors[-i]^2) / df);
        }
    }
    return(rstudentised);
}

#' @rdname plot.smooth
#' @export
plot.mes <- function(x, which=c(1,2,4,6), level=0.95, legend=FALSE,
                     ask=prod(par("mfcol")) < length(which) && dev.interactive(),
                     lowess=TRUE, ...){
    ellipsis <- list(...);

    # Define, whether to wait for the hit of "Enter"
    if(ask){
        oask <- devAskNewPage(TRUE);
        on.exit(devAskNewPage(oask));
    }

    # 1. Basic plot over time
    plot1 <- function(x, ...){
        yActuals <- actuals(x);
        if(!is.null(x$holdout)){
            yActuals <- ts(c(yActuals,x$holdout),start=start(yActuals),frequency=frequency(yActuals));
        }
        graphmaker(yActuals, x$forecast, fitted(x), main=x$model, legend=legend, parReset=FALSE, ...);
    }

    # 2 and 3: Standardised  / studentised residuals vs Fitted
    plot2 <- function(x, type="rstandard", ...){
        ellipsis <- list(...);

        ellipsis$x <- as.vector(fitted(x));
        if(type=="rstandard"){
            ellipsis$y <- as.vector(rstandard(x));
            yName <- "Standardised";
        }
        else{
            ellipsis$y <- as.vector(rstudent(x));
            yName <- "Studentised";
        }

        if(is.oes(x$occurrence)){
            ellipsis$x <- ellipsis$x[actuals(x$occurrence)!=0];
            ellipsis$y <- ellipsis$y[actuals(x$occurrence)!=0];
        }

        # Remove NAs
        if(any(is.na(ellipsis$x))){
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$x)];
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$y)];
        }

        if(!any(names(ellipsis)=="main")){
            ellipsis$main <- paste0(yName," Residuals vs Fitted");
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Fitted";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- paste0(yName," Residuals");
        }

        if(legend){
            if(ellipsis$x[length(ellipsis$x)]>mean(ellipsis$x)){
                legendPosition <- "bottomright";
            }
            else{
                legendPosition <- "topright";
            }
        }

        zValues <- switch(x$distribution,
                          "dlaplace"=qlaplace(c((1-level)/2, (1+level)/2), 0, 1),
                          "dalaplace"=qalaplace(c((1-level)/2, (1+level)/2), 0, 1, x$lambda),
                          "dlogis"=qlogis(c((1-level)/2, (1+level)/2), 0, 1),
                          "dt"=qt(c((1-level)/2, (1+level)/2), nobs(x)-nparam(x)),
                          "ds"=qs(c((1-level)/2, (1+level)/2), 0, 1),
                          # In the next one, the scale is debiased, taking n-k into account
                          "dinvgauss"=qinvgauss(c((1-level)/2, (1+level)/2), mean=1,
                                                dispersion=x$scale * nobs(x) / (nobs(x)-nparam(x))),
                          "dlnorm"=qlnorm(c((1-level)/2, (1+level)/2), 0, 1),
                          qnorm(c((1-level)/2, (1+level)/2), 0, 1));
        outliers <- which(ellipsis$y >zValues[2] | ellipsis$y <zValues[1]);
        # cat(paste0(round(length(outliers)/length(ellipsis$y),3)*100,"% of values are outside the bounds\n"));

        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- range(c(ellipsis$y,zValues), na.rm=TRUE);
            if(legend){
                if(legendPosition=="bottomright"){
                    ellipsis$ylim[1] <- ellipsis$ylim[1] - 0.2*diff(ellipsis$ylim);
                }
                else{
                    ellipsis$ylim[2] <- ellipsis$ylim[2] + 0.2*diff(ellipsis$ylim);
                }
            }
        }

        xRange <- range(ellipsis$x, na.rm=TRUE);
        xRange[1] <- xRange[1] - sd(ellipsis$x, na.rm=TRUE);
        xRange[2] <- xRange[2] + sd(ellipsis$x, na.rm=TRUE);

        do.call(plot,ellipsis);
        if(any(x$distribution==c("dlnorm","dinvgauss"))){
            abline(h=1, col="grey", lty=2);
        }
        else{
            abline(h=0, col="grey", lty=2);
        }
        polygon(c(xRange,rev(xRange)),c(zValues[1],zValues[1],zValues[2],zValues[2]),
                col="lightgrey", border=NA, density=10);
        abline(h=zValues, col="red", lty=2);
        if(length(outliers)>0){
            points(ellipsis$x[outliers], ellipsis$y[outliers], pch=16);
            text(ellipsis$x[outliers], ellipsis$y[outliers], labels=outliers, pos=4);
        }
        if(lowess){
            lines(lowess(ellipsis$x, ellipsis$y), col="red");
        }

        if(legend){
            if(lowess){
                legend(legendPosition,
                       legend=c(paste0(round(level,3)*100,"% bounds"),"outside the bounds","LOWESS line"),
                       col=c("red", "black","red"), lwd=c(1,NA,1), lty=c(2,1,1), pch=c(NA,16,NA));
            }
            else{
                legend(legendPosition,
                       legend=c(paste0(round(level,3)*100,"% bounds"),"outside the bounds"),
                       col=c("red", "black"), lwd=c(1,NA), lty=c(2,1), pch=c(NA,16));
            }
        }
    }

    # 4 and 5. Fitted vs |Residuals| or Fitted vs Residuals^2
    plot3 <- function(x, type="abs", ...){
        ellipsis <- list(...);

        ellipsis$x <- as.vector(fitted(x));
        if(type=="abs"){
            ellipsis$y <- abs(as.vector(residuals(x)));
        }
        else{
            ellipsis$y <- as.vector(residuals(x))^2;
        }

        if(is.oes(x$occurrence)){
            ellipsis$x <- ellipsis$x[ellipsis$y!=0];
            ellipsis$y <- ellipsis$y[ellipsis$y!=0];
        }
        # Remove NAs
        if(any(is.na(ellipsis$x))){
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$x)];
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$y)];
        }

        if(!any(names(ellipsis)=="main")){
            if(type=="abs"){
                ellipsis$main <- "|Residuals| vs Fitted";
            }
            else{
                ellipsis$main <- "Residuals^2 vs Fitted";
            }
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Fitted";
        }
        if(!any(names(ellipsis)=="ylab")){
            if(type=="abs"){
                ellipsis$ylab <- "|Residuals|";
            }
            else{
                ellipsis$ylab <- "Residuals^2";
            }
        }

        do.call(plot,ellipsis);
        if(any(x$distribution==c("dlnorm","dinvgauss"))){
            abline(h=1, col="grey", lty=2);
        }
        else{
            abline(h=0, col="grey", lty=2);
        }
        if(lowess){
            lines(lowess(ellipsis$x, ellipsis$y), col="red");
        }
    }

    # 6. Q-Q with the specified distribution
    plot4 <- function(x, ...){
        ellipsis <- list(...);

        ellipsis$y <- residuals(x);
        if(is.oes(x$occurrence)){
            ellipsis$y <- ellipsis$y[actuals(x$occurrence)!=0];
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Theoretical Quantile";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- "Actual Quantile";
        }

        if(any(x$distribution=="dnorm")){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ plot of normal distribution";
            }

            do.call(qqnorm, ellipsis);
            qqline(ellipsis$y);
        }
        else if(any(x$distribution=="dlnorm")){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ plot of log normal distribution";
            }
            ellipsis$x <- qlnorm(ppoints(500), meanlog=0, sdlog=x$scale);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qlnorm(p, meanlog=0, sdlog=x$scale));
        }
        else if(x$distribution=="dlaplace"){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of Laplace distribution";
            }
            ellipsis$x <- qlaplace(ppoints(500), mu=0, scale=x$scale);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qlaplace(p, mu=0, scale=x$scale));
        }
        else if(x$distribution=="dalaplace"){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- paste0("QQ-plot of Asymmetric Laplace with alpha=",round(x$lambda,3));
            }
            ellipsis$x <- qalaplace(ppoints(500), mu=0, scale=x$scale, alpha=x$lambda);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qalaplace(p, mu=0, scale=x$scale, alpha=x$lambda));
        }
        else if(x$distribution=="dlogis"){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of Logistic distribution";
            }
            ellipsis$x <- qlogis(ppoints(500), location=0, scale=x$scale);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qlogis(p, location=0, scale=x$scale));
        }
        else if(x$distribution=="ds"){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of S distribution";
            }
            ellipsis$x <- qs(ppoints(500), mu=0, scale=x$scale);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qs(p, mu=0, scale=x$scale));
        }
        else if(x$distribution=="dt"){
            # Standardise residuals
            ellipsis$y[] <- ellipsis$y / sd(ellipsis$y);
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of Student's distribution";
            }
            ellipsis$x <- qt(ppoints(500), df=x$scale);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qt(p, df=x$scale));
        }
        else if(x$distribution=="dinvgauss"){
            # Transform residuals for something meaningful
            # This is not 100% accurate, because the dispersion should change as well as mean...
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of Inverse Gaussian distribution";
            }
            ellipsis$x <- qinvgauss(ppoints(500), mean=1, dispersion=x$scale);

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) qinvgauss(p, mean=1, dispersion=x$scale));
        }
    }

    # 7 and 8. ACF and PACF
    plot5 <- function(x, type="acf", ...){
        ellipsis <- list(...);

        if(!any(names(ellipsis)=="main")){
            if(type=="acf"){
                ellipsis$main <- "Autocorrelation Function";
            }
            else{
                ellipsis$main <- "Partial Autocorrelation Function";
            }
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Lags";
        }
        if(!any(names(ellipsis)=="ylab")){
            if(type=="acf"){
                ellipsis$ylab <- "ACF";
            }
            else{
                ellipsis$ylab <- "PACF";
            }
        }

        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- c(-1,1);
        }

        if(type=="acf"){
            theValues <- acf(residuals(x), plot=FALSE, na.action=na.pass);
        }
        else{
            theValues <- pacf(residuals(x), plot=FALSE, na.action=na.pass);
        }
        ellipsis$x <- theValues$acf[-1];

        ellipsis$type <- "h"

        do.call(plot,ellipsis);
        abline(h=0, col="black", lty=1);
        abline(h=qnorm(c((1-level)/2, (1+level)/2),0,sqrt(1/nobs(x))), col="red", lty=2);
    }

    # 9. Plot of states
    plot6 <- function(x, ...){
        parDefault <- par(no.readonly = TRUE);
        if(any(unlist(gregexpr("C",x$model))==-1)){
            if(ncol(x$states)>10){
                message("Too many states. Plotting them one by one on several graphs.");
                if(is.null(ellipsis$main)){
                    ellipsisMain <- NULL;
                }
                else{
                    ellipsisMain <- ellipsis$main;
                }
                nPlots <- ceiling(ncol(x$states)/10);
                for(i in 1:nPlots){
                    if(is.null(ellipsisMain)){
                        ellipsis$main <- paste0("States of ",x$model,", part ",i);
                    }
                    ellipsis$x <- x$states[,(1+(i-1)*10):min(i*10,ncol(x$states))];
                    do.call(plot.ts, ellipsis);
                }
            }
            else{
                if(is.null(ellipsis$main)){
                    ellipsis$main <- paste0("States of ",x$model);
                }
                ellipsis$x <- x$states;
                do.call(plot.ts, ellipsis);
            }
        }
        else{
            # If we did combinations, we cannot return anything
            message("Combination of models was done. Sorry, but there is nothing to plot.");
        }
        par(parDefault);
    }

    # Do plots
    if(any(which==1)){
        plot1(x, ...);
    }

    if(any(which==2)){
        plot2(x, ...);
    }

    if(any(which==3)){
        plot2(x, "rstudent", ...);
    }

    if(any(which==4)){
        plot3(x, ...);
    }

    if(any(which==5)){
        plot3(x, type="squared", ...);
    }

    if(any(which==6)){
        plot4(x, ...);
    }

    if(any(which==7)){
        plot5(x, type="acf", ...);
    }

    if(any(which==8)){
        plot5(x, type="pacf", ...);
    }

    if(any(which==9)){
        plot6(x, ...);
    }
}

#' @export
pointLik.mes <- function(object, ...){
    distribution <- object$distribution;
    yInSample <- actuals(object);
    obsInSample <- nobs(object);
    if(is.oes(object$occurrence)){
        otLogical <- yInSample!=0;
        yFitted <- fitted(object) / fitted(object$occurrence);
    }
    else{
        otLogical <- rep(TRUE, obsInSample);
        yFitted <- fitted(object);
    }
    scale <- object$scale;

    likValues <- vector("numeric",obsInSample);
    likValues[otLogical] <- switch(distribution,
                                   "dnorm"=switch(Etype,
                                                  "A"=dnorm(x=yInSample[otLogical], mean=yFitted[otLogical],
                                                            sd=scale, log=TRUE),
                                                  "M"=dnorm(x=yInSample[otLogical], mean=yFitted[otLogical],
                                                            sd=scale*yFitted[otLogical], log=TRUE)),
                                   "dlogis"=switch(Etype,
                                                   "A"=dlogis(x=yInSample[otLogical], location=yFitted[otLogical],
                                                              scale=scale, log=TRUE),
                                                   "M"=dlogis(x=yInSample[otLogical], location=yFitted[otLogical],
                                                              scale=scale*yFitted[otLogical], log=TRUE)),
                                   "dlaplace"=switch(Etype,
                                                     "A"=dlaplace(q=yInSample[otLogical], mu=yFitted[otLogical],
                                                                  scale=scale, log=TRUE),
                                                     "M"=dlaplace(q=yInSample[otLogical], mu=yFitted[otLogical],
                                                                  scale=scale*yFitted[otLogical], log=TRUE)),
                                   "dt"=switch(Etype,
                                               "A"=dt(mesFitted$errors[otLogical], df=abs(lambda), log=TRUE),
                                               "M"=dt(mesFitted$errors[otLogical]*yFitted[otLogical],
                                                      df=abs(lambda), log=TRUE)),
                                   "ds"=switch(Etype,
                                               "A"=ds(q=yInSample[otLogical],mu=yFitted[otLogical],
                                                      scale=scale, log=TRUE),
                                               "M"=ds(q=yInSample[otLogical],mu=yFitted[otLogical],
                                                      scale=scale*sqrt(yFitted[otLogical]), log=TRUE)),
                                   "dalaplace"=switch(Etype,
                                                      "A"=dalaplace(q=yInSample[otLogical], mu=yFitted[otLogical],
                                                                    scale=scale, alpha=lambda, log=TRUE),
                                                      "M"=dalaplace(q=yInSample[otLogical], mu=yFitted[otLogical],
                                                                    scale=scale*yFitted[otLogical], alpha=lambda, log=TRUE)),
                                   "dlnorm"=dlnorm(x=yInSample[otLogical], meanlog=log(yFitted[otLogical]),
                                                   sdlog=scale, log=TRUE),
                                   "dinvgauss"=dinvgauss(x=yInSample[otLogical], mean=yFitted[otLogical],
                                                         dispersion=scale/yFitted[otLogical], log=TRUE))

    # If this is a mixture model, take the respective probabilities into account (differential entropy)
    if(is.oes(object$occurrence)){
        likValues[!otLogical] <- -switch(distribution,
                                         "dnorm" =,
                                         "dlnorm" = (log(sqrt(2*pi)*scale)+0.5),
                                         "dlogis" = 2,
                                         "dlaplace" =,
                                         "dalaplace" = (1 + log(2*scale)),
                                         "dt" = ((scale+1)/2 * (digamma((scale+1)/2)-digamma(scale/2)) +
                                                     log(sqrt(scale) * beta(scale/2,0.5))),
                                         "ds" = (2 + 2*log(2*scale)),
                                         "dinvgauss" = (0.5*(log(pi/2)+1+log(scale))));

        likValues[] <- likValues + pointLik(object$occurrence);
    }
    likValues <- ts(likValues, start=start(yFitted), frequency=frequency(yFitted));

    return(likValues);
}

#' @export
print.mes <- function(x, digits=4, ...){
    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"),2)," seconds"));
    cat(paste0("\nModel estimated: ",x$model));

    if(is.oes(x$occurrence)){
        occurrence <- switch(x$occurrence$occurrence,
                             "f"=,
                             "fixed"="Fixed probability",
                             "o"=,
                             "odds-ratio"="Odds ratio",
                             "i"=,
                             "inverse-odds-ratio"="Inverse odds ratio",
                             "d"=,
                             "direct"="Direct",
                             "g"=,
                             "general"="General");
        cat(paste0("\nOccurrence model type: ",occurrence));
    }

    distrib <- switch(x$distribution,
                      "dnorm" = "Normal",
                      "dlogis" = "Logistic",
                      "dlaplace" = "Laplace",
                      "dalaplace" = paste0("Asymmetric Laplace with lambda=",round(x$lambda,digits)),
                      "dt" = paste0("Student t with df=",round(x$lambda, digits)),
                      "ds" = "S",
                      "dlnorm" = "Log Normal",
                      # "dbcnorm" = paste0("Box-Cox Normal with lambda=",round(x$other$lambda,2)),
                      "dinvgauss" = "Inverse Gaussian"
    );
    if(is.oes(x$occurrence)){
        distrib <- paste0("Mixture of Bernoulli and ", distrib);
    }
    cat(paste0("\nDistribution assumed in the model: ", distrib));

    if(!is.null(x$persistence)){
        cat(paste0("\nPersistence vector g:\n"));
        if(is.matrix(x$persistence)){
            print(round(x$persistence,digits));
        }
        else{
            print(round(x$persistence,digits));
        }
    }

    if(!is.null(x$phi)){
        if(gregexpr("d",x$model)!=-1){
            cat(paste0("Damping parameter: ", round(x$phi,digits),"\n"));
        }
    }

    cat(paste0("\nLoss function type: ",x$loss));
    if(!is.null(x$lossValue)){
        cat(paste0("; Loss function value: ",round(x$lossValue,digits)));
        if(any(x$loss==c("LASSO","RIDGE"))){
            cat(paste0("; lambda=",x$lambda));
        }
    }

    cat("\nSample size: "); cat(nobs(x));
    cat("\nNumber of estimated parameters: "); cat(nparam(x));
    cat("\nNumber of degrees of freedom: "); cat(nobs(x)-nparam(x));

    if(x$loss=="likelihood" ||
       (any(x$loss==c("MSE","MSEh","MSCE")) & any(x$distribution==c("dnorm","dlnorm"))) ||
       (any(x$loss==c("MAE","MAEh","MACE")) & x$distribution=="dlaplace") ||
       (any(x$loss==c("HAM","HAMh","CHAM")) & x$distribution=="ds")){
           ICs <- c(AIC(x),AICc(x),BIC(x),BICc(x));
           names(ICs) <- c("AIC","AICc","BIC","BICc");
           cat("\nInformation criteria:\n");
           print(round(ICs,digits));
    }
    else{
        cat("\nInformation criteria are unavailable for the chosen loss & distribution.");
    }
}

#' @export
print.mesCombined <- function(x, digits=4, ...){
    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"),2)," seconds"));
    cat(paste0("\nModel estimated: ",x$model));

    cat(paste0("\nNumber of models combined: ", length(x$models)));

    cat("\nSample size: "); cat(nobs(x));
    cat("\nNumber of estimated parameters: "); cat(nparam(x));
    cat("\nNumber of degrees of freedom: "); cat(nobs(x)-nparam(x));
}

#' @export
print.mes.forecast <- function(x, ...){
    if(x$interval!="none"){
        returnedValue <- switch(x$side,
                                "both"=cbind(x$mean,x$lower,x$upper),
                                "lower"=cbind(x$mean,x$lower),
                                "upper"=cbind(x$mean,x$upper));
        colnames(returnedValue) <- switch(x$side,
                                          "both"=c("Point forecast",
                                                   paste0("Lower bound (",mean((1-x$level)/2)*100,"%)"),
                                                   paste0("Upper bound (",mean((1+x$level)/2)*100,"%)")),
                                          "lower"=c("Point forecast",
                                                   paste0("Lower bound (",mean((1-x$level))*100,"%)")),
                                          "upper"=c("Point forecast",
                                                   paste0("Upper bound (",mean(x$level)*100,"%)")));
    }
    else{
        returnedValue <- x$mean;
    }
    print(returnedValue);
}

#' @importFrom stats sigma
#' @export
sigma.mes <- function(object, ...){
    return(sqrt(switch(object$distribution,
                       "dnorm"=,
                       "dlogis"=,
                       "dlaplace"=,
                       "dt"=,
                       "ds"=,
                       "dalaplace"=sum(residuals(object)^2),
                       "dlnorm"=sum(log(residuals(object))^2),
                       "dinvgauss"=sum((residuals(object)-1)^2))
                /(nobs(object)-nparam(object))));
}

#' @export
summary.mes <- function(object, level=0.95, ...){
    ourReturn <- list(timeElapsed=object$timeElapsed, model=object$model);

    occurrence <- NULL;
    if(is.oes(object$occurrence)){
        occurrence <- switch(object$occurrence$occurrence,
                             "f"=,
                             "fixed"="Fixed probability",
                             "o"=,
                             "odds-ratio"="Odds ratio",
                             "i"=,
                             "inverse-odds-ratio"="Inverse odds ratio",
                             "d"=,
                             "direct"="Direct",
                             "g"=,
                             "general"="General");
    }
    ourReturn$occurrence <- occurrence;
    ourReturn$distribution <- object$distribution;

    # Collect parameters and their standard errors
    parametersConfint <- confint(object, level=level);
    parametersConfint[,2:3] <- coef(object) + parametersConfint[,2:3];
    parametersTable <- cbind(coef(object),parametersConfint);
    rownames(parametersTable) <- names(coef(object));
    colnames(parametersTable) <- c("Estimate","Std. Error",
                                   paste0("Lower ",(1-level)/2*100,"%"),
                                   paste0("Upper ",(1+level)/2*100,"%"));
    ourReturn$coefficients <- parametersTable;
    ourReturn$loss <- object$loss;
    ourReturn$lossValue <- object$lossValue;
    ourReturn$nobs <- nobs(object);
    ourReturn$nparam <- nparam(object);

    if(object$loss=="likelihood" ||
       (any(object$loss==c("MSE","MSEh","MSCE")) & any(object$distribution==c("dnorm","dlnorm"))) ||
       (any(object$loss==c("MAE","MAEh","MACE")) & object$distribution=="dlaplace") ||
       (any(object$loss==c("HAM","HAMh","CHAM")) & object$distribution=="ds")){
        ICs <- c(AIC(object),AICc(object),BIC(object),BICc(object));
        names(ICs) <- c("AIC","AICc","BIC","BICc");
        ourReturn$ICs <- ICs;
    }
    return(structure(ourReturn, class="summary.mes"));
}

print.summary.mes <- function(x, ...){
    ellipsis <- list(...);
    if(!any(names(ellipsis)=="digits")){
        digits <- 4;
    }
    else{
        digits <- ellipsis$digits;
    }

    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"), digits)," seconds"));
    cat(paste0("\nModel estimated: ",x$model));

    if(!is.null(x$occurrence)){
        cat(paste0("\nOccurrence model type: ",x$occurrence));
    }

    distrib <- switch(x$distribution,
                      "dnorm" = "Normal",
                      "dlogis" = "Logistic",
                      "dlaplace" = "Laplace",
                      "dalaplace" = paste0("Asymmetric Laplace with lambda=",round(x$lambda,digits)),
                      "dt" = paste0("Student t with df=",round(x$lambda, digits)),
                      "ds" = "S",
                      "dlnorm" = "Log Normal",
                      # "dbcnorm" = paste0("Box-Cox Normal with lambda=",round(x$other$lambda,2)),
                      "dinvgauss" = "Inverse Gaussian"
    );
    if(!is.null(x$occurrence)){
        distrib <- paste0("Mixture of Bernoulli and ", distrib);
    }
    cat(paste0("\nDistribution assumed in the model: ", distrib));

    if(!is.null(x$coefficients)){
        cat("\nCoefficients:\n");
        print(round(x$coefficients,digits));
    }

    cat(paste0("\nLoss function type: ",x$loss));
    if(!is.null(x$lossValue)){
        cat(paste0("; Loss function value: ",round(x$lossValue,digits)));
        if(any(x$loss==c("LASSO","RIDGE"))){
            cat(paste0("; lambda=",x$lambda));
        }
    }

    cat("\nSample size: "); cat(x$nobs);
    cat("\nNumber of estimated parameters: "); cat(x$nparam);
    cat("\nNumber of degrees of freedom: "); cat(x$nobs-x$nparam);

    if(x$loss=="likelihood" ||
       (any(x$loss==c("MSE","MSEh","MSCE")) & any(x$distribution==c("dnorm","dlnorm"))) ||
       (any(x$loss==c("MAE","MAEh","MACE")) & x$distribution=="dlaplace") ||
       (any(x$loss==c("HAM","HAMh","CHAM")) & x$distribution=="ds")){
        cat("\n");
        print(round(x$ICs,digits));
    }
    else{
        cat("\nInformation criteria are unavailable for the chosen loss & distribution.\n");
    }
}

# Work in progress...
# predict.mes <- function(object, newxreg=NULL, interval=c("none", "confidence", "prediction"),
#                         level=0.95, side=c("both","upper","lower"), ...){
#
#     h <- nrow(newxreg);
#     lagsModelAll <- object$lags;
#     componentsNumber <- length(lagsModelAll);
#     componentsNumberSeasonal <- sum(lagsModelAll>1);
#     lagsModelMax <- max(lagsModelAll);
#
#     matVt <- t(object$states[obsStates-(lagsModelMax:1)+1,,drop=FALSE]);
#     matWt <- tail(object$measurement,h);
#     matF <- object$transition;
#     vecG <- object$persistence;
#
#     model <- modelType(test);
#     Etype <- errorType(object);
#     Ttype <- substr(model,2,2);
#     Stype <- substr(model,nchar(model),nchar(model));
#
#     mesForecast <- mesForecasterWrap(matVt, tail(matWt,h), matF, vecG,
#                                      lagsModelAll, Etype, Ttype, Stype,
#                                      componentsNumber, componentsNumberSeasonal, h);
#
#     yForecast <- mesForecast$yForecast;
#
#     return(yForecast);

# Work in progress...
#' @param nsim Number of iterations to do in case of \code{interval="simulated"}.
#' @rdname forecast.smooth
#' @importFrom stats rnorm rlogis rt rlnorm qnorm qlogis qt qlnorm
#' @importFrom statmod rinvgauss qinvgauss
#' @importFrom greybox rlaplace rs ralaplace qlaplace qs qalaplace
#' @export
forecast.mes <- function(object, h=10, newxreg=NULL,
                         interval=c("none", "simulated", "approximate", "semiparametric", "nonparametric"),
                         level=0.95, side=c("both","upper","lower"), cumulative=FALSE, nsim=10000, ...){

    ellipsis <- list(...);

    # If the horizon is zero, just construct fitted and potentially confidence interval thingy
    if(h<=0){
        return(predict(object, newxreg=newxreg,
                       interval=interval,
                       level=level, side=side, ...))
    }
    interval <- match.arg(interval[1],c("none", "simulated", "approximate", "semiparametric", "nonparametric","parametric"));
    if(interval=="parametric"){
        warning("The parameter 'interval' does not accept 'parametric' anymore. We use 'approximate' value instead.",
                call.=FALSE)
        interval <- "approximate";
    }
    side <- match.arg(side);

    # Technical parameters
    lagsModelAll <- lagsModel <- lags(object);
    componentsNumber <- length(lagsModel);
    componentsNumberSeasonal <- sum(lagsModel>1);
    lagsModelMax <- max(lagsModel);
    obsStates <- nrow(object$states);
    obsInSample <- nobs(object);

    # Model type
    model <- modelType(object);
    Etype <- errorType(object);
    Ttype <- substr(model,2,2);
    Stype <- substr(model,nchar(model),nchar(model));

    # ts structure
    yForecastStart <- time(actuals(object))[obsInSample]+deltat(actuals(object));
    yFrequency <- frequency(actuals(object));

    # All the important matrices
    matVt <- t(object$states[obsStates-(lagsModelMax:1)+1,,drop=FALSE]);
    matWt <- tail(object$measurement,h);
    vecG <- matrix(object$persistence, ncol=1);
    if(!is.null(object$xreg)){
        xregNumber <- ncol(object$xreg);
        lagsModelAll <- rbind(lagsModelAll,rep(1,xregNumber));
        matWt[,componentsNumber+c(1:xregNumber)] <- newxreg[1:h,];
        vecG <- rbind(vecG,object$xregPersistence);
    }
    else{
        xregNumber <- 0;
    }
    matF <- diag(componentsNumber+xregNumber);
    matF[1:componentsNumber,1:componentsNumber] <- object$transition;

    # Produce point forecasts
    mesForecast <- mesForecasterWrap(matVt, matWt, matF, vecG,
                                     lagsModelAll, Etype, Ttype, Stype,
                                     componentsNumber, componentsNumberSeasonal, h);

    # If this is a mixture model, produce forecasts for the occurrence
    if(!is.null(object$occurrence)){
        occurrenceModel <- TRUE;
        pForecast <- forecast(object$occurrence,h=h)$mean;
    }
    else{
        occurrenceModel <- FALSE;
        pForecast <- rep(1, h);
    }

    # Cumulative forecasts have only one observation
    if(cumulative){
        yForecast <- yUpper <- yLower <- ts(vector("numeric", 1), start=yForecastStart, frequency=yFrequency);
        yForecast[] <- sum(mesForecast$yForecast * pForecast);
    }
    else{
        yForecast <- yUpper <- yLower <- ts(vector("numeric", h), start=yForecastStart, frequency=yFrequency);
        yForecast[] <- mesForecast$yForecast * pForecast;
    }

    if(interval!="none"){
        # If this is an occurrence model, then take probability into account in the level.
        if(occurrenceModel){
            levelNew <- (level-(1-pForecast))/pForecast;
            levelNew[levelNew<0] <- 0;
        }
        else{
            levelNew <- level;
        }

        if(cumulative){
            levelLow <- levelUp <- vector("numeric",1);
        }
        else{
            levelLow <- levelUp <- vector("numeric",h);
        }
        if(side=="both"){
            levelLow[] <- (1-levelNew)/2;
            levelUp[] <- (1+levelNew)/2;
        }
        else if(side=="upper"){
            levelLow[] <- rep(0,length(levelNew));
            levelUp[] <- levelNew;
        }
        else{
            levelLow[] <- 1-levelNew;
            levelUp[] <- rep(1,length(levelNew));
        }
        levelLow[levelLow<0] <- 0;
        levelUp[levelUp<0] <- 0;
    }

    # If simulated intervals are needed...
    if(interval=="simulated"){
        arrVt <- array(NA, c(componentsNumber+xregNumber, h+lagsModelMax, nsim));
        arrVt[,1:lagsModelMax,] <- rep(matVt,nsim);
        sigmaValue <- sigma(object);
        matErrors <- matrix(switch(object$distribution,
                                   "dnorm"=rnorm(h*nsim, 0, sigmaValue),
                                   "dlogis"=rlogis(h*nsim, 0, sigmaValue*sqrt(3)/pi),
                                   "dlaplace"=rlaplace(h*nsim, 0, sigmaValue/2),
                                   "dt"=rt(h*nsim, nobs(object)-nparam(object)),
                                   "ds"=rs(h*nsim, 0, (sigmaValue^2/120)^0.25),
                                   "dalaplace"=ralaplace(h*nsim, 0,
                                                         sqrt(sigmaValue^2*object$lambda^2*(1-object$lambda)^2/(object$lambda^2+(1-object$lambda)^2)),
                                                         object$lambda),
                                   "dlnorm"=rlnorm(h*nsim, 0, sigmaValue)-1,
                                   "dinvgauss"=rinvgauss(h*nsim, 1, dispersion=sigmaValue^2)-1,
                                   ),
                            h,nsim);
        # This stuff is needed in order to produce adequate values for weird models
        EtypeModified <- Etype;
        if(Etype=="A" && any(object$distribution==c("dlnorm","dinvgauss"))){
            EtypeModified[] <- "M";
        }
        matOt <- matrix(rbinom(h*nsim, 1, pForecast), h, nsim);
        matG <- matrix(vecG, componentsNumber+xregNumber, nsim);

        ySimulated <- mesSimulatorwrap(arrVt, matErrors, matOt, array(matF,c(dim(matF),nsim)), matWt, matG,
                                       EtypeModified, Ttype, Stype, lagsModelAll,
                                       componentsNumberSeasonal, componentsNumber)$matrixYt;

        #### Note that the cumulative doesn't work with oes at the moment!
        if(cumulative){
            yForecast[] <- mean(colSums(ySimulated,na.rm=T)*pForecast);
            yLower[] <- quantile(colSums(ySimulated,na.rm=T),levelLow,type=7);
            yUpper[] <- quantile(colSums(ySimulated,na.rm=T),levelUp,type=7);
        }
        else{
            # yForecast[] <- apply(ySimulated,1,mean,na.rm=T) * pForecast;
            for(i in 1:h){
                yLower[i] <- quantile(ySimulated[i,],levelLow[i],na.rm=T,type=7);
                yUpper[i] <- quantile(ySimulated[i,],levelUp[i],na.rm=T,type=7);
            }
        }
    }
    # This option will use h-steps ahead variance and produce intervals for it
    # This will rely on multicov() method
    else if(interval=="approximate"){
        s2 <- sigma(object)^2;
        vcovMulti <- covarAnal(lagsModelAll, h, matWt[1,,drop=FALSE], matF, vecG, s2);

        # Do either the variance of sum, or a diagonal
        if(cumulative){
            vcovMulti <- sum(vcovMulti);
        }
        else{
            vcovMulti <- diag(vcovMulti);
        }

        if(object$distribution=="dnorm"){
            if(errorType(object)=="A"){
                yLower[] <- qnorm(levelLow, yForecast, sqrt(vcovMulti));
                yUpper[] <- qnorm(levelUp, yForecast, sqrt(vcovMulti));
            }
            else{
                yLower[] <- yForecast*qnorm(levelLow, 1, sqrt(vcovMulti));
                yUpper[] <- yForecast*qnorm(levelUp, 1, sqrt(vcovMulti));
            }
        }
        if(object$distribution=="dlogis"){
            if(errorType(object)=="A"){
                yLower[] <- qlogis(levelLow, yForecast, sqrt(vcovMulti*3)/pi);
                yUpper[] <- qlogis(levelUp, yForecast, sqrt(vcovMulti*3)/pi);
            }
            else{
                yLower[] <- yForecast*qlogis(levelLow, 1, sqrt(vcovMulti*3)/pi);
                yUpper[] <- yForecast*qlogis(levelUp, 1, sqrt(vcovMulti*3)/pi);
            }
        }
        if(object$distribution=="dlaplace"){
            if(errorType(object)=="A"){
                yLower[] <- qlaplace(levelLow, yForecast, sqrt(vcovMulti/2));
                yUpper[] <- qlaplace(levelUp, yForecast, sqrt(vcovMulti/2));
            }
            else{
                yLower[] <- yForecast*qlaplace(levelLow, 1, sqrt(vcovMulti/2));
                yUpper[] <- yForecast*qlaplace(levelUp, 1, sqrt(vcovMulti/2));
            }
        }
        if(object$distribution=="dt"){
            df <- nobs(object) - nparam(object);
            if(errorType(object)=="A"){
                yLower[] <- yForecast + sqrt(vcovMulti)*qt(levelLow, df);
                yUpper[] <- yForecast + sqrt(vcovMulti)*qt(levelUp, df);
            }
            else{
                yLower[] <- yForecast*(1 + sqrt(vcovMulti)*qt(levelLow, df));
                yUpper[] <- yForecast*(1 + sqrt(vcovMulti)*qt(levelUp, df));
            }
        }
        if(object$distribution=="ds"){
            if(errorType(object)=="A"){
                yLower[] <- qs(levelLow, yForecast, (vcovMulti/120)^0.25);
                yUpper[] <- qs(levelUp, yForecast, (vcovMulti/120)^0.25);
            }
            else{
                yLower[] <- yForecast*qs(levelLow, 1, (vcovMulti/120)^0.25);
                yUpper[] <- yForecast*qs(levelUp, 1, (vcovMulti/120)^0.25);
            }
        }
        if(object$distribution=="dalaplace"){
            lambda <- object$lambda;
            if(errorType(object)=="A"){
                yLower[] <- qalaplace(levelLow, yForecast,
                                      sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
                yUpper[] <- qalaplace(levelUp, yForecast,
                                      sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
            }
            else{
                yLower[] <- yForecast*qalaplace(levelLow, 1,
                                                sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
                yUpper[] <- yForecast*qalaplace(levelUp, 1,
                                                sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
            }
        }
        if(object$distribution=="dlnorm"){
            yLower[] <- yForecast*qlnorm(levelLow, 0, sqrt(vcovMulti));
            yUpper[] <- yForecast*qlnorm(levelUp, 0, sqrt(vcovMulti));
        }
        else if(object$distribution=="dinvgauss"){
            if(errorType(object)=="A"){
                yLower[] <- yForecast*qinvgauss(levelLow, 1, dispersion=vcovMulti);
                yUpper[] <- yForecast*qinvgauss(levelUp, 1, dispersion=vcovMulti);
            }
            else{
                vcovMulti <- mesVarAnal(h, matWt[1,,drop=FALSE], vecG, s2);
                yLower[] <- yForecast*qinvgauss(levelLow, 1, dispersion=vcovMulti);
                yUpper[] <- yForecast*qinvgauss(levelUp, 1, dispersion=vcovMulti);
            }
        }
    }
    # This option will extract the matrix of multisteps errors from mes and build intervals based on that
    # This will rely on rmultistep() method
    # else if(any(interval==c("semiparametric","nonparametric"))){
    # }
    else{
        yUpper[] <- yLower[] <- NA;
    }

    # Fix of prediction intervals depending on what has happened
    if(interval!="none"){
        # Make sensible values out of those weird quantiles
        if(!cumulative){
            if(Etype=="A"){
                yLower[levelLow==0] <- -Inf;
            }
            else{
                yLower[levelLow==0] <- 0;
            }
            yUpper[levelUp==1] <- Inf;
        }
        else{
            if(Etype=="A" && (levelLow==0)){
                yLower[] <- -Inf;
            }
            else if(Etype=="M" && (levelLow==0)){
                yLower[] <- 0;
            }
            if(levelUp==1){
                yUpper[] <- Inf;
            }
        }

        # Check what we have from the occurrence model
        if(occurrenceModel){
            # If there are NAs, then there's no variability and no intervals.
            if(any(is.na(yUpper))){
                yUpper[is.na(yUpper)] <- (yForecast/pForecast)[is.na(yUpper)];
            }
            if(any(is.na(yLower))){
                yLower[is.na(yLower)] <- 0;
            }
        }
    }

    return(structure(list(mean=yForecast, lower=yLower, upper=yUpper, model=object,
                          level=level, interval=interval, side=side, cumulative=cumulative),
                     class=c("mes.forecast","smooth.forecast","forecast")));
}

forecast.mesCombined <- function(object, h=10, newxreg=NULL,
                                 interval=c("none", "simulated", "approximate", "semiparametric", "nonparametric"),
                                 level=0.95, side=c("both","upper","lower"), cumulative=FALSE, nsim=10000, ...){
    interval <- match.arg(interval);
    side <- match.arg(side);

    # ts structure
    yForecastStart <- time(actuals(object))[nobs(object)]+deltat(actuals(object));
    yFrequency <- frequency(actuals(object));

    # Cumulative forecasts have only one observation
    if(cumulative){
        yForecast <- yUpper <- yLower <- ts(vector("numeric", 1), start=yForecastStart, frequency=yFrequency);
    }
    else{
        yForecast <- yUpper <- yLower <- ts(vector("numeric", h), start=yForecastStart, frequency=yFrequency);
    }

    # The list contains 8 elements
    mesForecasts <- vector("list",8);
    names(mesForecasts)[c(1:3)] <- c("mean","lower","upper");
    for(i in 1:length(object$models)){
        mesForecasts[] <- forecast.mes(object$models[[i]], h=h, newxreg=newxreg,
                                     interval=interval,
                                     level=level, side=side, cumulative=cumulative, nsim=nsim, ...);
        yForecast[] <- yForecast + mesForecasts$mean * object$ICw[i];
        yUpper[] <- yUpper + mesForecasts$upper * object$ICw[i];
        yLower[] <- yLower + mesForecasts$lower * object$ICw[i];
    }

    # Get rid of specific models
    object$models <- NULL;

    return(structure(list(mean=yForecast, lower=yLower, upper=yUpper, model=object,
                          level=level, interval=interval, side=side, cumulative=cumulative),
                     class=c("mes.forecast","smooth.forecast","forecast")));
}

#' @export
multicov.mes <- function(object, type=c("analytical","empirical","simulated"), ...){
    h <- length(object$holdout);
    lagsModelAll <- lags(object);
    componentsNumber <- length(lagsModelAll);
    s2 <- sigma(object)^2;
    matWt <- tail(object$measurement,h);
    vecG <- matrix(object$persistence, ncol=1);
    if(!is.null(object$xreg)){
        xregNumber <- ncol(object$xreg);
        lagsModelAll <- rbind(lagsModelAll,rep(1,xregNumber));
        matWt[,componentsNumber+c(1:xregNumber)] <- newxreg[1:h,];
        vecG <- rbind(vecG,object$xregPersistence);
    }
    else{
        xregNumber <- 0;
    }
    matF <- diag(componentsNumber+xregNumber);
    matF[1:componentsNumber,1:componentsNumber] <- object$transition;

    covarMat <- covarAnal(lagsModelAll, h, matWt[1,,drop=FALSE], matF, vecG, s2);

    return(covarMat);
}

#' @export
vcov.mes <- function(object, ...){
    modelReturn <- mes(actuals(object), model=object, FI=TRUE);
    return(solve(modelReturn$FI));
}

##### Other functions to implement #####
# rmultistep.mes <- function(object, ...){}
# accuracy.mes <- function(object, holdout, ...){}
# simulate.mes <- function(object, nsim=1, seed=NULL, obs=NULL, ...){}
