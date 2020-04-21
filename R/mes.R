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
#' \item dllaplace - Log Laplace distribution,
#' \item ds - Log S distribution,
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
#' MES(M,Ad,M)[m1,m2,...], where m1, m2, ... are the lags specified by the
#' \code{lags} parameter.
#' There are several options for the \code{model} besides the conventional ones,
#' which rely on information criteria:
#' \enumerate{
#' \item \code{model="ZZZ"} means that the model will be selected based on the
#' chosen information criteria type. The Branch and Bound is used in the process.
#' \item \code{model="XXX"} means that only additive components are tested, using
#' Branch and Bound.
#' \item \code{model="YYY"} implies selecting between multiplicative components.
#' \item \code{model="CCC"} trigers the combination of forecasts of models using
#' information criteria weights (Kolassa, 2011).
#' \item combinations between these four and the classical components are also
#' accepted. For example, \code{model="CAY"} will combine models with additive
#' trend and either none or multiplicative seasonality.
#' \item \code{model="PPP"} will produce the selection between pure additive and
#' pure multiplicative models. "P" stands for "Pure". This cannot be mixed with
#' other types of components.
#' \item \code{model="FFF"} will select between all the 30 types of models. "F"
#' stands for "Full". This cannot be mixed with other types of components.
#' \item The parameter \code{model} can also be a vector of names of models for a
#' finer tuning (pool of models). For example, \code{model=c("ANN","AAA")} will
#' estimate only two models and select the best of them.
#' }
#'
#' Also, \code{model} can accept a previously estimated mes model and use all
#' its parameters.
#'
#' Keep in mind that model selection with "Z" components uses Branch and Bound
#' algorithm and may skip some models that could have slightly smaller
#' information criteria. If you want to do a exhaustive search, you would need
#' to list all the models to check as a vector.
#'
#' The default value is set to \code{"ZXZ"}, because the multiplicative trend is explosive
#' and dangerous. It should be used only for each separate time series, not for the
#' atomated predictions for big  datasets.
#'
#' @param lags Defines lags for the corresponding components. All components
#' count, starting from level, so ETS(M,M,M) model for monthly data will have
#' lags=c(1,1,12). If fractional numbers are provided, then it is assumed that
#' the data is not periodic. The parameter \code{date} is then needed in order
#' to setup the appropriate time series structure.
#' @param orders The order of ARIMA to be included in the model. This should be passed
#' either as a vector (in which case the non-seasonal ARIMA is assumed) or as a list of
#' a type \code{orders=list(ar=c(p,P),i=c(d,D),ma=c(q,Q))}, in which case the \code{lags}
#' variable is used in order to determine the seasonality m. See \link[smooth]{msarima}
#' for details.
#' @param formula Formula to use in case of explanatory variables. If \code{NULL},
#' then all the variables are used as is. Only considered if \code{xreg} is not
#' \code{NULL} and \code{xregDo="use"}.
#' @param distribution what density function to assume for the error term. The full
#' name of the distribution should be provided, starting with the letter "d" -
#' "density". The names align with the names of distribution functions in R.
#' For example, see \link[stats]{dnorm}. For detailed explanation of available
#' distributions, see vignette in greybox package: \code{vignette("greybox","alm")}.
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
#' are available: \code{MAEh}, \code{TMAE}, \code{GTMAE}, \code{MACE},
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
#' Also, a model produced using \link[smooth]{oes} or \link[greybox]{alm} function
#' can be used here.
#' @param ic The information criterion to use in the model selection / combination
#' procedure.
#' @param bounds The type of bounds for the persistence to use in the model
#' estimation. Can be either \code{admissible} - guaranteeing the stability of the
#' model, \code{traditional} - restricting the values with (0, 1) or \code{none} - no
#' restrictions (potentially dangerous).
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
#' @param silent Specifies, whether to provide the progress of the function or not.
#' If \code{TRUE}, then the function will print what it does and how much it has
#' already done.
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
#' \item \code{maxeval} - maximum number of evaluations to carry out (default is 200 for
#' data with lags <= 12 and 1000 for the larger lags);
#' \item \code{maxtime} - stop, when the optimisation time (in seconds) exceeds this;
#' \item \code{xtol_rel} - the relative precision of the optimiser (the default is 1E-6);
#' \item \code{xtol_abs} - the absolute precision of the optimiser (the default is 0 -
#' not used);
#' \item \code{algorithm} - the algorithm to use in optimisation
#' (\code{"NLOPT_LN_BOBYQA"} by default);
#' \item \code{print_level} - the level of output for the optimiser (0 by default).
#' }
#' You can read more about these parameters by running the function
#' \link[nloptr]{nloptr.print.options}.
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
#' \item \code{formula} - the formula used for the explanatory variables expansion,
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
#' @importFrom forecast forecast na.interp
#' @importFrom greybox dlaplace dalaplace ds stepwise alm is.occurrence
#' @importFrom stats dnorm dlogis dt dlnorm frequency
#' @importFrom statmod dinvgauss
#' @importFrom nloptr nloptr
#' @importFrom pracma hessian
#' @useDynLib mes
#' @export mes
mes <- function(y, model="ZXZ", lags=c(frequency(y)), orders=list(ar=c(0),i=c(0),ma=c(0)), formula=NULL,
                distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                               "dlnorm","dllaplace","dls","dinvgauss"),
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
    # paste0() is needed in order to get rid of potential issues with names
    responseName <- paste0(deparse(substitute(y)),collapse="");

    #### Check the parameters of the function and create variables based on them ####
    parametersChecker(y, model, lags, formula, orders,
                      persistence, phi, initial,
                      distribution, loss, h, holdout, occurrence, ic, bounds,
                      xreg, xregDo, xregInitial, xregPersistence, responseName,
                      silent, modelDo, ParentEnvironment=environment(), ellipsis, fast=FALSE);
    # Remove xreg if it was provided, just to preserve some memory
    rm(xreg);

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
                        obsStates, obsInSample, obsAll, componentsNumber, componentsNumberSeasonal,
                        componentsNames, otLogical,
                        yInSample, persistence, persistenceEstimate, phi,
                        initialValue, initialType,
                        xregExist, xregInitialsProvided, xregPersistence,
                        xregModel, xregData, xregNumber, xregNames){
        # Matrix of states. Time in columns, components in rows
        matVt <- matrix(NA, componentsNumber+xregNumber, obsStates, dimnames=list(c(componentsNames,xregNames),NULL));

        # Measurement rowvector
        matWt <- matrix(1, obsAll, componentsNumber+xregNumber, dimnames=list(NULL,c(componentsNames,xregNames)));
        # If xreg are provided, then fill in the respective values in Wt vector
        if(xregExist){
            matWt[,componentsNumber+1:xregNumber] <- xregData;
        }

        # Transition matrix
        matF <- diag(componentsNumber+xregNumber);

        # Persistence vector
        vecG <- matrix(0, componentsNumber+xregNumber, 1, dimnames=list(c(componentsNames,xregNames),NULL));
        if(!persistenceEstimate){
            vecG[1:componentsNumber,] <- persistence;
        }
        if(xregExist){
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
                # If either Etype or Stype are multiplicative, do multiplicative decomposition
                decompositionType <- c("additive","multiplicative")[any(c(Etype,Stype)=="M")+1];
                yDecomposition <- msdecompose(yInSample, lags[lags!=1], type=decompositionType);
                j <- 1;
                # level
                matVt[j,1:lagsModelMax] <- mean(yInSample[1:lagsModelMax]);
                j <- j+1;
                if(Ttype!="N"){
                    if(Ttype=="A" && Stype=="M"){
                        # level fix
                        matVt[j-1,1:lagsModelMax] <- exp(mean(log(yInSample[otLogical][1:lagsModelMax])));
                        # trend
                        matVt[j,1:lagsModelMax] <- prod(yDecomposition$initial)-yDecomposition$initial[1];
                    }
                    else if(Ttype=="M" && Stype=="A"){
                        # level fix
                        matVt[j-1,1:lagsModelMax] <- exp(mean(log(yInSample[otLogical][1:lagsModelMax])));
                        # trend
                        matVt[j,1:lagsModelMax] <- sum(yDecomposition$initial)/yDecomposition$initial[1];
                    }
                    else{
                        # trend
                        matVt[j,1:lagsModelMax] <- yDecomposition$initial[2];
                    }
                    # This is a failsafe for multiplicative trend models, so that the thing does not explode
                    if(Ttype=="M" && matVt[j,1:lagsModelMax]>1.1){
                        matVt[j,1:lagsModelMax] <- 1;
                    }
                    j <- j+1;
                }
                #### Seasonal components
                # For pure models use stuff as is
                if(all(c(Etype,Stype)=="A") || all(c(Etype,Stype)=="M") ||
                   (Etype=="A" & Stype=="M")){
                    for(i in 1:componentsNumberSeasonal){
                        matVt[i+j-1,(lagsModelMax-lagsModel[i+j-1])+1:lagsModel[i+j-1]] <- yDecomposition$seasonal[[i]];
                    }
                }
                # For mixed models use a different set of initials
                else if(Etype=="M" && Stype=="A"){
                    for(i in 1:componentsNumberSeasonal){
                        matVt[i+j-1,(lagsModelMax-lagsModel[i+j-1])+
                                  1:lagsModel[i+j-1]] <- log(yDecomposition$seasonal[[i]])*min(yInSample[otLogical]);
                    }
                }
            }
            # Non-seasonal models
            else{
                matVt[1,1] <- mean(yInSample[1:max(lagsModelMax,ceiling(obsInSample*0.2))]);
                if(Ttype!="N"){
                    matVt[2,1] <- switch(Ttype,
                                         "A" = mean(diff(yInSample[1:max(lagsModelMax+1,ceiling(obsInSample*0.2))]),na.rm=TRUE),
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
            if(Etype=="A" || xregInitialsProvided || is.null(xregModel[[2]])){
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
                # A special type of model which is not safe: AAM, MAA, MAM
                if((Etype=="A" && Ttype=="A" && Stype=="M") || (Etype=="A" && Ttype=="M" && Stype=="A") ||
                   ((initialType=="backcasting") &&
                    ((Etype=="M" && Ttype=="A" && Stype=="A") || (Etype=="M" && Ttype=="A" && Stype=="M")))){
                    B[j:componentsNumber] <- c(0.01,0,rep(0,componentsNumberSeasonal))[j:componentsNumber];
                }
                # MMA is the worst. Set everything to zero and see if anything can be done...
                else if((Etype=="M" && Ttype=="M" && Stype=="A")){
                    B[j:componentsNumber] <- c(0,0,rep(0,componentsNumberSeasonal))[j:componentsNumber];
                }
                else if(Etype=="M" && Ttype=="A"){
                    if(initialType=="backcasting"){
                        B[j:componentsNumber] <- c(0.1,0,rep(0.11,componentsNumberSeasonal))[j:componentsNumber];
                    }
                    else{
                        B[j:componentsNumber] <- c(0.1,0.05,rep(0.11,componentsNumberSeasonal))[j:componentsNumber];
                    }
                }
                else{
                    if(Ttype!="N"){
                        B[j:componentsNumber] <- c(0.1,0.05,rep(0.11,componentsNumberSeasonal))[j:componentsNumber];
                    }
                    else{
                        B[j:componentsNumber] <- c(0.1,rep(0.2,componentsNumberSeasonal))[j:componentsNumber];
                    }
                }
            }
            else{
                if(Ttype!="N"){
                    B[j:componentsNumber] <- c(0.1,0.05,rep(0.11,componentsNumberSeasonal))[j:componentsNumber];
                }
                else{
                    B[j:componentsNumber] <- c(0.1,rep(0.2,componentsNumberSeasonal))[j:componentsNumber];
                }
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
                if(componentsNumberSeasonal>1){
                    for(k in i:componentsNumber){
                        B[j+0:(lagsModel[k]-1)] <- matVt[k,(lagsModelMax-lagsModel[k])+1:lagsModel[k]];
                        names(B)[j+0:(lagsModel[k]-1)] <- paste0("seasonal",k-i+1,"_",1:lagsModel[k]);
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
                else{
                    B[j+0:(lagsModel[componentsNumber]-1)] <- matVt[componentsNumber,(lagsModelMax-lagsModel[componentsNumber])+
                                                                        1:lagsModel[componentsNumber]];
                    names(B)[j+0:(lagsModel[componentsNumber]-1)] <- paste0("seasonal_",1:lagsModel[componentsNumber]);
                    if(Stype=="A"){
                        Bl[j+0:(lagsModel[componentsNumber]-1)] <- -Inf;
                        Bu[j+0:(lagsModel[componentsNumber]-1)] <- Inf;
                    }
                    else{
                        Bl[j+0:(lagsModel[componentsNumber]-1)] <- 0;
                        Bu[j+0:(lagsModel[componentsNumber]-1)] <- Inf;
                    }
                    j <- j+lagsModel[componentsNumber];
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
                        "dllaplace"=switch(Etype,
                                        "A"=sum(abs(log(1+errors/yFitted))/obsInSample),
                                        "M"=sum(abs(log(1+errors))/obsInSample)),
                        "dls"=switch(Etype,
                                        "A"=sum(sqrt(abs(log(1+errors/yFitted)))/obsInSample),
                                        "M"=sum(sqrt(abs(log(1+errors)))/obsInSample)),
                        "dinvgauss"=switch(Etype,
                                           "A"=sum((errors/yFitted)^2/(1+errors/yFitted))/obsInSample,
                                           "M"=sum((errors)^2/(1+errors))/obsInSample),
                                           # "M"=mean((errors)^2/(1+errors))),
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
                   xregExist, xregInitialsEstimate, xregPersistenceEstimate,
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

        # Check the bounds, classical restrictions
        if(bounds=="usual"){
            if(any(mesElements$vecG>1) || any(mesElements$vecG<0)){
                return(1E+300);
            }
            if(Ttype!="N"){
                if((mesElements$vecG[2]>mesElements$vecG[1])){
                    return(1E+300);
                }
                if(Stype!="N" && any(mesElements$vecG[-c(1,2)]>(1-mesElements$vecG[1]))){
                    return(1E+300);
                }
            }
            else{
                if(Stype!="N" && any(mesElements$vecG[-1]>(1-mesElements$vecG[1]))){
                    return(1E+300);
                }
            }
            # This is the restriction on the damping parameter
            if(phiEstimate && (mesElements$matF[2,2]>1 || mesElements$matF[2,2]<0)){
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
                                       "dllaplace"=dlaplace(q=log(yInSample[otLogical]), mu=log(mesFitted$yFitted[otLogical]),
                                                            scale=scale, log=TRUE)-log(yInSample[otLogical]),
                                       "dls"=ds(q=log(yInSample[otLogical]), mu=log(mesFitted$yFitted[otLogical]),
                                                scale=scale, log=TRUE)-log(yInSample[otLogical]),
                                       # "dinvgauss"=dinvgauss(x=1+mesFitted$errors, mean=1,
                                       #                       dispersion=scale, log=TRUE)));
                                       # "dinvgauss"=dinvgauss(x=yInSampleNew, mean=mesFitted$yFitted,
                                       #                       dispersion=scale/mesFitted$yFitted, log=TRUE)));
                                       "dinvgauss"=dinvgauss(x=yInSample[otLogical], mean=mesFitted$yFitted[otLogical],
                                                             dispersion=scale/mesFitted$yFitted[otLogical], log=TRUE)));

                # Differential entropy for the logLik of occurrence model
                if(occurrenceModel || any(!otLogical)){
                    CFValue <- CFValue + switch(distribution,
                                                "dnorm" =,
                                                "dlnorm" = obsZero*(log(sqrt(2*pi)*scale)+0.5),
                                                "dlogis" = obsZero*2,
                                                "dlaplace" =,
                                                "dllaplace" =,
                                                "dalaplace" = obsZero*(1 + log(2*scale)),
                                                "dt" = obsZero*((scale+1)/2 *
                                                                    (digamma((scale+1)/2)-digamma(scale/2)) +
                                                                    log(sqrt(scale) * beta(scale/2,0.5))),
                                                "ds" =,
                                                "dls" = obsZero*(2 + 2*log(2*scale)),
                                                # "dinvgauss" = obsZero*(0.5*(log(pi/2)+1+suppressWarnings(log(scale)))));
                                                # "dinvgauss" =0);
                                                "dinvgauss" = 0.5*(obsZero*(log(pi/2)+1+suppressWarnings(log(scale)))-
                                                                       sum(log(mesFitted$yFitted[!otLogical]))));
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
            mesErrors <- mesErrorerWrap(mesFitted$matVt, mesElements$matWt, mesElements$matF,
                                        lagsModelAll, Etype, Ttype, Stype,
                                        componentsNumber, componentsNumberSeasonal, h,
                                        yInSample, ot);

            # This is a fix for the multistep in case of Etype=="M", assuming logN
            if(Etype=="M"){
                mesErrors[] <- log(1+mesErrors);
            }

            # Not done yet: "aMSEh","aTMSE","aGTMSE","aMSCE","aGPL"
            CFValue <- switch(loss,
                              "MSEh"=sum(mesErrors[,h]^2)/(obsInSample-h),
                              "TMSE"=sum(colSums(mesErrors^2)/(obsInSample-h)),
                              "GTMSE"=sum(log(colSums(mesErrors^2)/(obsInSample-h))),
                              "MSCE"=sum(rowSums(mesErrors)^2)/(obsInSample-h),
                              "MAEh"=sum(abs(mesErrors[,h]))/(obsInSample-h),
                              "TMAE"=sum(colSums(abs(mesErrors))/(obsInSample-h)),
                              "GTMAE"=sum(log(colSums(abs(mesErrors))/(obsInSample-h))),
                              "MACE"=sum(abs(rowSums(mesErrors)))/(obsInSample-h),
                              "HAMh"=sum(sqrt(abs(mesErrors[,h])))/(obsInSample-h),
                              "THAM"=sum(colSums(sqrt(abs(mesErrors)))/(obsInSample-h)),
                              "GTHAM"=sum(log(colSums(sqrt(abs(mesErrors)))/(obsInSample-h))),
                              "CHAM"=sum(sqrt(abs(rowSums(mesErrors))))/(obsInSample-h),
                              "GPL"=log(det(t(mesErrors) %*% mesErrors/(obsInSample-h))),
                              0);

        }

        if(is.na(CFValue) || is.nan(CFValue)){
            CFValue[] <- 1e+300;
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
                          xregExist, xregInitialsEstimate, xregPersistenceEstimate,
                          xregNumber,
                          bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate){

        if(!multisteps){
            if(any(loss==c("LASSO","RIDGE"))){
                return(0);
            }
            else{
                distributionNew <- switch(loss,
                                          "MSE"=switch(Etype,"A"="dnorm","M"="dlnorm"),
                                          "MAE"=switch(Etype,"A"="dlaplace","M"="dllaplace"),
                                          "HAM"=switch(Etype,"A"="ds","M"="dls"),
                                          distribution);
                logLikReturn <- -CF(B,
                                    Etype, Ttype, Stype, yInSample,
                                    ot, otLogical, occurrenceModel, obsInSample,
                                    componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                                    matVt, matWt, matF, vecG, componentsNumberSeasonal,
                                    persistenceEstimate, phiEstimate, initialType,
                                    xregExist, xregInitialsEstimate, xregPersistenceEstimate,
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
            logLikReturn <- CF(B,
                               Etype, Ttype, Stype, yInSample,
                               ot, otLogical, occurrenceModel, obsInSample,
                               componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                               matVt, matWt, matF, vecG, componentsNumberSeasonal,
                               persistenceEstimate, phiEstimate, initialType,
                               xregExist, xregInitialsEstimate, xregPersistenceEstimate,
                               xregNumber,
                               bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);

            logLikReturn[] <- -switch(loss,
                                      "MSEh"=, "aMSEh"=, "TMSE"=, "aTMSE"=, "MSCE"=, "aMSCE"=
                                          (obsInSample-h)/2*(log(2*pi)+1+log(logLikReturn)),
                                      "GTMSE"=, "aGTMSE"=
                                          (obsInSample-h)/2*(log(2*pi)+1+logLikReturn),
                                      "MAEh"=, "TMAE"=, "GTMAE"=, "MACE"=
                                          (obsInSample-h)*(log(2)+1+log(logLikReturn)),
                                      "HAMh"=, "THAM"=, "GTHAM"=, "CHAM"=
                                          (obsInSample-h)*(log(4)+2+2*log(logLikReturn)),
                                      #### Divide GPL by 8 in order to make it comparable with the univariate ones
                                      "GPL"=, "aGPL"=
                                          (obsInSample-h)/2*(h*log(2*pi)+h+logLikReturn)/h);

            # This is not well motivated at the moment, but should make likelihood comparable, taking T instead of T-h
            logLikReturn[] <- logLikReturn / (obsInSample-h) * obsInSample;

            # In case of multiplicative model, we assume log- distribution
            if(Etype=="M"){
                logLikReturn[] <- logLikReturn - sum(log(yInSample));
            }

            return(logLikReturn);
        }
    }

    #### The function estimates the ETS model and returns B, logLik, nParam and CF(B) ####
    estimator <- function(Etype, Ttype, Stype, lags,
                          obsStates, obsInSample,
                          yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                          initialType, initialValue,
                          xregExist, xregInitialsProvided, xregInitialsEstimate,
                          xregPersistence, xregPersistenceEstimate,
                          xregModel, xregData, xregNumber, xregNames, xregDo,
                          ot, otLogical, occurrenceModel, pFitted,
                          bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate){

        # Create the basic variables
        mesArchitect <- architector(Etype, Ttype, Stype, lags, xregNumber, obsInSample, initialType);
        list2env(mesArchitect, environment());

        # Create the matrices for the specific ETS model
        mesCreated <- creator(Etype, Ttype, Stype,
                              lags, lagsModel, lagsModelMax, lagsLength, lagsModelAll,
                              obsStates, obsInSample, obsAll, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregExist, xregInitialsProvided, xregPersistence,
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
        if(distribution=="default"){
            distributionNew <- switch(Etype,
                                     "A"=switch(loss,
                                                "MAEh"=, "MACE"=, "MAE"="dlaplace",
                                                "HAMh"=, "CHAM"=, "HAM"="ds",
                                                "MSEh"=, "MSCE"=, "MSE"=, "GPL"=, "likelihood"=, "dnorm"),
                                     "M"=switch(loss,
                                                "MAEh"=, "MACE"=, "MAE"="dllaplace",
                                                "HAMh"=, "CHAM"=, "HAM"="dls",
                                                "MSEh"=, "MSCE"=, "MSE"=, "GPL"="dlnorm",
                                                "likelihood"=, "dinvgauss"));
        }
        else{
            distributionNew <- distribution;
        }
        # print(B)
        # print(Etype)
        # print(Ttype)
        # print(Stype)

        # Parameters are chosen to speed up the optimisation process and have decent accuracy
        res <- suppressWarnings(nloptr(B, CF, lb=lb, ub=ub,
                                       opts=list(algorithm=algorithm, xtol_rel=xtol_rel, xtol_abs=xtol_abs,
                                                 maxeval=maxeval, maxtime=maxtime, print_level=print_level),
                                       Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                                       ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, obsInSample=obsInSample,
                                       componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                                       matVt=mesCreated$matVt, matWt=mesCreated$matWt, matF=mesCreated$matF, vecG=mesCreated$vecG,
                                       componentsNumberSeasonal=componentsNumberSeasonal,
                                       persistenceEstimate=persistenceEstimate, phiEstimate=phiEstimate, initialType=initialType,
                                       xregExist=xregExist, xregInitialsEstimate=xregInitialsEstimate,
                                       xregPersistenceEstimate=xregPersistenceEstimate, xregNumber=xregNumber,
                                       bounds=bounds, loss=loss, distribution=distributionNew, horizon=horizon, multisteps=multisteps,
                                       lambda=lambda, lambdaEstimate=lambdaEstimate));

        if(is.infinite(res$objective) || res$objective==1e+300){
            # If the optimisation didn't work, give it another try with zero initials for smoothing parameters
            B[1:componentsNumber] <- 0;
            res <- suppressWarnings(nloptr(B, CF, lb=lb, ub=ub,
                                           opts=list(algorithm=algorithm, xtol_rel=xtol_rel, maxeval=maxeval,
                                                     maxtime=maxtime, print_level=print_level),
                                           Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                                           ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, obsInSample=obsInSample,
                                           componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                                           matVt=mesCreated$matVt, matWt=mesCreated$matWt, matF=mesCreated$matF, vecG=mesCreated$vecG,
                                           componentsNumberSeasonal=componentsNumberSeasonal,
                                           persistenceEstimate=persistenceEstimate, phiEstimate=phiEstimate, initialType=initialType,
                                           xregExist=xregExist, xregInitialsEstimate=xregInitialsEstimate,
                                           xregPersistenceEstimate=xregPersistenceEstimate, xregNumber=xregNumber,
                                           bounds=bounds, loss=loss, distribution=distributionNew, horizon=horizon, multisteps=multisteps,
                                           lambda=lambda, lambdaEstimate=lambdaEstimate));
        }

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
                                              xregExist, xregInitialsEstimate, xregPersistenceEstimate,
                                              xregNumber,
                                              bounds, loss, distributionNew, horizon, multisteps, lambda, lambdaEstimate),
                                    nobs=obsInSample,df=nParamEstimated,class="logLik");

        #### If we do variables selection, do it here, then reestimate the model. ####
        if(xregDo=="select"){
            # Fill in the matrices
            mesElements <- filler(B,
                                  Ttype, Stype, componentsNumber, lagsModel, lagsModelMax,
                                  mesCreated$matVt, mesCreated$matWt, mesCreated$matF, mesCreated$vecG,
                                  persistenceEstimate, phiEstimate, initialType,
                                  xregInitialsEstimate, xregPersistenceEstimate, xregNumber);

            # Fit the model to the data
            mesFitted <- mesFitterWrap(mesElements$matVt, mesElements$matWt, mesElements$matF, mesElements$vecG,
                                       lagsModelAll, Etype, Ttype, Stype, componentsNumber, componentsNumberSeasonal,
                                       yInSample, ot, initialType=="backcasting");

            # Extract the errors and amend them to correspond to the distribution
            errors <- mesFitted$errors+switch(Etype,"A"=0,"M"=1);

            # Call the xregSelector providing the original matrix with the data
            if(Etype=="A"){
                xregModel[[1]] <- xregSelector(errors=errors, xregData=xregDataOriginal, ic=ic,
                                               df=length(B)+1, distribution=distributionNew, occurrence=oesModel);
                xregNumber <- length(xregModel[[1]]$xregInitial);
                xregNames <- names(xregModel[[1]]$xregInitial);
            }
            else{
                xregModel[[2]] <- xregSelector(errors=errors, xregData=xregDataOriginal, ic=ic,
                                               df=length(B)+1, distribution=distributionNew, occurrence=oesModel);
                xregNumber <- length(xregModel[[2]]$xregInitial);
                xregNames <- names(xregModel[[2]]$xregInitial);
            }

            # If there are some variables, then do the proper reestimation and return the new values
            if(xregNumber>0){
                xregExist[] <- TRUE;
                xregInitialsEstimate[] <- TRUE;
                xregPersistenceEstimate[] <- TRUE;
                xregData <- xregDataOriginal[,xregNames,drop=FALSE];

                return(estimator(Etype, Ttype, Stype, lags,
                                 obsStates, obsInSample,
                                 yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                 initialType, initialValue,
                                 xregExist, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames, xregDo="use",
                                 ot, otLogical, occurrenceModel, pFitted,
                                 bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate));

            }
        }

        return(list(B=B, CFValue=CFValue, nParamEstimated=nParamEstimated, logLikMESValue=logLikMESValue,
                    xregExist=xregExist, xregData=xregData, xregNumber=xregNumber, xregNames=xregNames, xregModel=xregModel,
                    xregInitialsEstimate=xregInitialsEstimate, xregPersistenceEstimate=xregPersistenceEstimate));
    }


    #### The function creates a pool of models and selects the best of them ####
    selector <- function(model, modelsPool, allowMultiplicative,
                         Etype, Ttype, Stype, damped, lags,
                         obsStates, obsInSample,
                         yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                         initialType, initialValue,
                         xregExist, xregInitialsProvided, xregInitialsEstimate,
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
            # Align error and seasonality, if the error was not forced to be additive
            if(any(substr(poolSmall,3,3)=="M") && all(Etype!=c("A","X"))){
                multiplicativeSeason <- (substr(poolSmall,3,3)=="M");
                poolSmall[multiplicativeSeason] <- paste0("M",substr(poolSmall[multiplicativeSeason],2,3));
            }
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
                    phi[] <- 1;
                    phiEstimate[] <- FALSE;
                    Stype[] <- substring(modelCurrent,3,3);
                }

                results[[i]] <- estimator(Etype, Ttype, Stype, lags,
                                          obsStates, obsInSample,
                                          yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                          initialType, initialValue,
                                          xregExist, xregInitialsProvided, xregInitialsEstimate,
                                          xregPersistence, xregPersistenceEstimate,
                                          xregModel, xregData, xregNumber, xregNames, xregDo,
                                          ot, otLogical, occurrenceModel, pFitted,
                                          bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);
                results[[i]]$IC <- ICFunction(results[[i]]$logLikMESValue);
                results[[i]]$Etype <- Etype;
                results[[i]]$Ttype <- Ttype;
                results[[i]]$Stype <- Stype;
                results[[i]]$phiEstimate <- phiEstimate;
                results[[j]]$phi <- phi;
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
            # print(modelCurrent)
            Etype <- substring(modelCurrent,1,1);
            Ttype <- substring(modelCurrent,2,2);
            if(nchar(modelCurrent)==4){
                phi[] <- 0.95;
                Stype <- substring(modelCurrent,4,4);
                phiEstimate <- TRUE;
            }
            else{
                phi[] <- 1;
                Stype <- substring(modelCurrent,3,3);
                phiEstimate <- FALSE;
            }

            results[[j]] <- estimator(Etype, Ttype, Stype, lags,
                                      obsStates, obsInSample,
                                      yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                      initialType, initialValue,
                                      xregExist, xregInitialsProvided, xregInitialsEstimate,
                                      xregPersistence, xregPersistenceEstimate,
                                      xregModel, xregData, xregNumber, xregNames, xregDo,
                                      ot, otLogical, occurrenceModel, pFitted,
                                      bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);
            results[[j]]$IC <- ICFunction(results[[j]]$logLikMESValue);
            results[[j]]$Etype <- Etype;
            results[[j]]$Ttype <- Ttype;
            results[[j]]$Stype <- Stype;
            results[[j]]$phiEstimate <- phiEstimate;
            if(phiEstimate){
                results[[j]]$phi <- results[[j]]$B[names(results[[j]]$B)=="phi"];
            }
            else{
                results[[j]]$phi <- 1;
            }
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

    ##### Function uses residuals in order to determine the needed xreg #####
    xregSelector <- function(errors, xregData, ic, df, distribution, occurrence){
        stepwiseModel <- stepwise(cbind(as.data.frame(errors),xregData[1:obsInSample,,drop=FALSE]), ic=ic, df=df,
                                  distribution=distribution, occurrence=occurrence, silent=TRUE);
        return(list(xregInitial=coef(stepwiseModel)[-1],other=stepwiseModel$other));
    }

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

        # Produce forecasts if the horizon is non-zero
        if(horizon>0){
            yForecast <- ts(rep(NA, horizon), start=yForecastStart, frequency=dataFreq);
            yForecast[] <- mesForecasterWrap(matVt[,obsInSample+(1:lagsModelMax),drop=FALSE], tail(matWt,horizon), matF,
                                             lagsModelAll, Etype, Ttype, Stype,
                                             componentsNumber, componentsNumberSeasonal, horizon);
            #### Make safety checks
            # If there are NaN values
            if(any(is.nan(yForecast))){
                yForecast[is.nan(yForecast)] <- 0;
            }
            # If there are negative values in the multiplicative model
            # if((Etype=="M") && any(yForecast<=0)){
            #     yForecast[yForecast<=0] <- 0.01;
            # }

            # Amend forecasts, multiplying by probability
            if(occurrenceModel && !occurrenceModelProvided){
                yForecast[] <- yForecast * c(forecast(oesModel, h=h)$mean);
            }
            else if(occurrenceModel && occurrenceModelProvided){
                yForecast[] <- yForecast * pForecast;
            }
        }
        else{
            yForecast <- ts(NA, start=yForecastStart, frequency=dataFreq);
        }

        # If the distribution is default, change it according to the error term
        if(distribution=="default"){
            distribution[] <- switch(Etype,
                                     "A"=switch(loss,
                                                "MAEh"=, "MACE"=, "MAE"="dlaplace",
                                                "HAMh"=, "CHAM"=, "HAM"="ds",
                                                "MSEh"=, "MSCE"=, "GPL"=, "MSE"=,
                                                "aMSEh"=, "aMSCE"=, "aGPL"=, "likelihood"=, "dnorm"),
                                     "M"="dinvgauss");
            if(multisteps && Etype=="M"){
                distribution[] <- switch(loss,
                                         "MAEh"=, "MACE"=, "MAE"="dllaplace",
                                         "HAMh"=, "CHAM"=, "HAM"="dls",
                                         "MSEh"=, "MSCE"=, "GPL"=, "MSE"=,
                                         "aMSEh"=, "aMSCE"=, "aGPL"=, "dlnorm");
            }
        }

        if(initialType=="optimal"){
            initialValue <- vector("numeric", sum(lagsModel));
            j <- 0;
            for(i in 1:length(lagsModel)){
                initialValue[j+1:lagsModel[i]] <- tail(matVt[i,1:lagsModelMax],lagsModel[i]);
                j <- j + lagsModel[i];
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

        if(xregInitialsEstimate){
            xregInitial <- matVt[-c(1:componentsNumber),lagsModelMax];
        }

        scale <- scaler(distribution, Etype, errors[otLogical], yFitted[otLogical], obsInSample, lambda);
        yFitted <- ts(yFitted, start=dataStart, frequency=dataFreq);

        return(list(model=NA, timeElapsed=NA,
                    y=NA, holdout=NA, fitted=yFitted, residuals=ts(errors, start=dataStart, frequency=dataFreq),
                    forecast=yForecast, states=ts(t(matVt), start=(time(y)[1] - deltat(y)*lagsModelMax),
                                                  frequency=dataFreq),
                    persistence=persistence, phi=phi, transition=matF,
                    measurement=matWt, initialType=initialType, initial=initialValue,
                    nParam=parametersNumber, occurrence=oesModel, xreg=xregData,
                    xregInitial=xregInitial, xregPersistence=xregPersistence, formula=formula,
                    loss=loss, lossValue=CFValue, logLik=logLikMESValue, distribution=distribution,
                    scale=scale, lambda=lambda, B=B, lags=lagsModel, FI=FI));
    }

    #### Deal with occurrence model ####
    if(occurrenceModel && !occurrenceModelProvided){
        oesModel <- suppressWarnings(oes(ot, model=model, occurrence=occurrence, ic=ic, h=horizon,
                                         holdout=FALSE, bounds="usual", xreg=xregData, xregDo=xregDo, silent=TRUE));
        pFitted[] <- fitted(oesModel);
        parametersNumber[1,3] <- nparam(oesModel);
        # print(oesModel)
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

    ##### Prepare stuff for the variables selection if xregDo="select" #####
    if(xregDo=="select"){
        # First, record the original parameters
        xregExistOriginal <- xregExist;
        xregInitialsProvidedOriginal <- xregInitialsProvided;
        xregInitialsEstimateOriginal <- xregInitialsEstimate;
        xregPersistenceOriginal <- xregPersistence;
        xregPersistenceProvidedOriginal <- xregPersistenceProvided;
        xregPersistenceEstimateOriginal <- xregPersistenceEstimate;
        xregModelOriginal <- xregModel;
        xregDataOriginal <- xregData;
        xregNumberOriginal <- xregNumber;
        xregNamesOriginal <- xregNames;

        # Set the parameters to zero and do simple ETS
        xregExist[] <- FALSE;
        xregInitialsProvided <- FALSE;
        xregInitialsEstimate[] <- FALSE;
        xregPersistence <- 0;
        xregPersistenceProvided <- FALSE;
        xregPersistenceEstimate[] <- FALSE;
        xregModel[[1]] <- NULL;
        xregModel[[2]] <- NULL;
        xregData <- NULL;
        xregNumber[] <- 0;
        xregNames <- NULL;
    }

    ##### Estimate the specified model #####
    if(modelDo=="estimate"){
        # Estimate the parameters of the demand sizes model
        mesEstimated <- estimator(Etype, Ttype, Stype, lags,
                                 obsStates, obsInSample,
                                 yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                 initialType, initialValue,
                                 xregExist, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames, xregDo,
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
                              obsStates, obsInSample, obsAll, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregExist, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);
        list2env(mesCreated, environment());

        icSelection <- ICFunction(mesEstimated$logLikMESValue);

        ####!!! If the occurrence is auto, then compare this with the model with no occurrence !!!####

        parametersNumber[1,1] <- (sum(lagsModel)*(initialType=="optimal") + phiEstimate +
                                      componentsNumber*persistenceEstimate + xregNumber*xregInitialsEstimate +
                                      xregNumber*xregPersistenceEstimate + 1);
        if(xregExist){
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
                                 xregExist, xregInitialsProvided, xregInitialsEstimate,
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
                              obsStates, obsInSample, obsAll, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregExist, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);
        list2env(mesCreated, environment());

        parametersNumber[1,1] <- (sum(lagsModel)*(initialType=="optimal") + phiEstimate +
                                      componentsNumber*persistenceEstimate + xregNumber*xregInitialsEstimate +
                                      xregNumber*xregPersistenceEstimate + 1);
        if(xregExist){
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
                                 xregExist, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames,
                                 ot, otLogical, occurrenceModel, pFitted, ICFunction,
                                 bounds, loss, distribution, horizon, multisteps, lambda, lambdaEstimate);

        icSelection <- mesSelected$icSelection;

        icBest <- min(icSelection);
        mesSelected$icWeights  <- (exp(-0.5*(icSelection-icBest)) /
                                       sum(exp(-0.5*(icSelection-icBest))));

        # This is a failsafe mechanism, just to make sure that the ridiculous models don't impact forecasts
        mesSelected$icWeights[mesSelected$icWeights<1e-5] <- 0
        mesSelected$icWeights <- mesSelected$icWeights/sum(mesSelected$icWeights);

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
                                  obsStates, obsInSample, obsAll, componentsNumber, componentsNumberSeasonal,
                                  componentsNames, otLogical,
                                  yInSample, persistence, persistenceEstimate, phi,
                                  initialValue, initialType,
                                  xregExist, xregInitialsProvided, xregPersistence,
                                  xregModel, xregData, xregNumber, xregNames);

            mesSelected$results[[i]]$matVt <- mesCreated$matVt;
            mesSelected$results[[i]]$matWt <- mesCreated$matWt;
            mesSelected$results[[i]]$matF <- mesCreated$matF;
            mesSelected$results[[i]]$vecG <- mesCreated$vecG;

            parametersNumber[1,1] <- (sum(lagsModel)*(initialType=="optimal") + phiEstimate +
                                          componentsNumber*persistenceEstimate + xregNumber*xregInitialsEstimate +
                                          xregNumber*xregPersistenceEstimate + 1);
            if(xregExist){
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
        if(distribution=="default"){
            distributionNew <- switch(Etype,
                                      "A"=switch(loss,
                                                 "MAEh"=, "MACE"=, "MAE"="dlaplace",
                                                 "HAMh"=, "CHAM"=, "HAM"="ds",
                                                 "MSEh"=, "MSCE"=, "MSE"=, "GPL"=, "likelihood"=, "dnorm"),
                                      "M"=switch(loss,
                                                 "MAEh"=, "MACE"=, "MAE"="dllaplace",
                                                 "HAMh"=, "CHAM"=, "HAM"="dls",
                                                 "MSEh"=, "MSCE"=, "MSE"=, "GPL"="dlnorm",
                                                 "likelihood"=, "dinvgauss"));
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
                              obsStates, obsInSample, obsAll, componentsNumber, componentsNumberSeasonal,
                              componentsNames, otLogical,
                              yInSample, persistence, persistenceEstimate, phi,
                              initialValue, initialType,
                              xregExist, xregInitialsProvided, xregPersistence,
                              xregModel, xregData, xregNumber, xregNames);
        list2env(mesCreated, environment());

        CFValue <- CF(B=0, Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                      ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, obsInSample=obsInSample,
                      componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                      matVt=mesCreated$matVt, matWt=mesCreated$matWt, matF=mesCreated$matF, vecG=mesCreated$vecG,
                      componentsNumberSeasonal=componentsNumberSeasonal,
                      persistenceEstimate=persistenceEstimate, phiEstimate=phiEstimate, initialType=initialType,
                      xregExist=xregExist, xregInitialsEstimate=xregInitialsEstimate,
                      xregPersistenceEstimate=xregPersistenceEstimate, xregNumber=xregNumber,
                      bounds=bounds, loss=loss, distribution=distributionNew, horizon=horizon, multisteps=multisteps,
                      lambda=lambda, lambdaEstimate=lambdaEstimate);

        parametersNumber[1,1] <- parametersNumber[1,4] <- 1;
        logLikMESValue <- structure(logLikMES(B=0,
                                              Etype, Ttype, Stype, yInSample,
                                              ot, otLogical, occurrenceModel, pFitted, obsInSample,
                                              componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                                              mesCreated$matVt, mesCreated$matWt, mesCreated$matF, mesCreated$vecG, componentsNumberSeasonal,
                                              persistenceEstimate, phiEstimate, initialType,
                                              xregExist, xregInitialsEstimate, xregPersistenceEstimate,
                                              xregNumber,
                                              bounds, loss, distributionNew, horizon, multisteps, lambda, lambdaEstimate)
                                    ,nobs=obsInSample,df=parametersNumber[1,4],class="logLik")

        icSelection <- ICFunction(logLikMESValue);
        # If Fisher Information is required, do that analytically
        if(FI){
            # If B is not provided, then use the standard thing
            if(is.null(B)){
                BValues <- initialiser(Etype, Ttype, Stype, componentsNumberSeasonal,
                                       componentsNumber, lagsModel, lagsModelMax, mesCreated$matVt,
                                       TRUE, damped, "optimal",
                                       xregExist, FALSE, xregNumber, FALSE);
                # Create the vector of initials for the optimisation
                B <- BValues$B;
            }

            # Define parameters just for FI calculation
            if(initialType=="provided" && any(names(B)=="level")){
                initialTypeFI <- "optimal";
            }
            else{
                initialTypeFI <- initialType;
            }

            #### This will not work in cases, when initials for xreg were provided from the beginning!
            if(xregExist && initialTypeFI=="optimal"){
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

            FI <- hessian(logLikMES, B, Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                          ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel, pFitted=pFitted, obsInSample=obsInSample,
                          componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                          matVt=matVt, matWt=matWt, matF=matF, vecG=vecG,
                          componentsNumberSeasonal=componentsNumberSeasonal,
                          persistenceEstimate=persistenceEstimateFI, phiEstimate=phiEstimateFI, initialType=initialTypeFI,
                          xregExist=xregExist, xregInitialsEstimate=xregInitialsEstimateFI,
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
    yInSample <- ts(yInSample,start=dataStart, frequency=dataFreq);
    if(holdout){
        yHoldout <- ts(yHoldout, start=yForecastStart, frequency=dataFreq);
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
        if(any(yNAValues)){
            modelReturned$y[yNAValues[1:obsInSample]] <- NA;
            if(length(yNAValues)==obsAll){
                modelReturned$holdout[yNAValues[-c(1:obsInSample)]] <- NA;
            }
            modelReturned$residuals[yNAValues[1:obsInSample]] <- NA;
        }

        class(modelReturned) <- c("mes","smooth");
    }
    #### Return the combined model ####
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
            modelReturned$models[[i]]$fitted[is.na(modelReturned$models[[i]]$fitted)] <- 0;
            yFittedCombined[] <- yFittedCombined + modelReturned$models[[i]]$fitted * mesSelected$icWeights[i];
            if(h>0){
                modelReturned$models[[i]]$forecast[is.na(modelReturned$models[[i]]$forecast)] <- 0;
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
            parametersNumberOverall[1,1] <- parametersNumber[1,1] + parametersNumber[1,1] * mesSelected$icWeights[i];
            modelReturned$models[[i]]$y <- yInSample;
            if(any(yNAValues)){
                modelReturned$models[[i]]$y[yNAValues[1:obsInSample]] <- NA;
                if(length(yNAValues)==obsAll){
                    modelReturned$models[[i]]$holdout[yNAValues[-c(1:obsInSample)]] <- NA;
                }
                modelReturned$models[[i]]$residuals[yNAValues[1:obsInSample]] <- NA;
            }

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
        modelReturned$fitted <- ts(yFittedCombined,start=dataStart, frequency=dataFreq);
        modelReturned$residuals <- yInSample - yFittedCombined;
        if(any(yNAValues)){
            modelReturned$y[yNAValues[1:obsInSample]] <- NA;
            if(length(yNAValues)==obsAll){
                modelReturned$holdout[yNAValues[-c(1:obsInSample)]] <- NA;
            }
            modelReturned$residuals[yNAValues[1:obsInSample]] <- NA;
        }
        modelReturned$forecast <- ts(yForecastCombined,start=yForecastStart, frequency=dataFreq);
        parametersNumberOverall[1,4] <- sum(parametersNumberOverall[1,1:3]);
        modelReturned$nParam <- parametersNumberOverall;
        modelReturned$ICw <- mesSelected$icWeights;
        # These two are needed just to make basic methods work
        modelReturned$distribution <- distribution;
        modelReturned$scale <- sqrt(mean(modelReturned$residuals^2,na.rm=TRUE));
        class(modelReturned) <- c("mesCombined","mes","smooth");
    }
    modelReturned$ICs <- icSelection;

    # Error measures if there is a holdout
    if(holdout){
        modelReturned$accuracy <- measures(yHoldout,modelReturned$forecast,yInSample);
    }

    if(!silent){
        plot(modelReturned, 7);
    }

    return(modelReturned);
}

#### Technical methods ####
#' @export
lags.mes <- function(object, ...){
    if(!is.null(object$xreg)){
        return(c(object$lags,rep(1,ncol(object$xreg))));
    }
    else{
        return(object$lags);
    }
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

    # 1. Fitted vs Actuals values
    plot1 <- function(x, ...){
        ellipsis <- list(...);

        # Get the actuals and the fitted values
        ellipsis$y <- c(actuals(x));
        if(is.occurrence(x)){
            if(any(x$distribution==c("plogis","pnorm"))){
                ellipsis$y <- (ellipsis$y!=0)*1;
            }
        }
        ellipsis$x <- c(fitted(x));

        # If this is a mixture model, remove zeroes
        if(is.occurrence(x$occurrence)){
            ellipsis$x <- ellipsis$x[ellipsis$y!=0];
            ellipsis$y <- ellipsis$y[ellipsis$y!=0];
        }

        # Remove NAs
        if(any(is.na(ellipsis$x))){
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$x)];
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$x)];
        }
        if(any(is.na(ellipsis$y))){
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$y)];
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$y)];
        }

        # Title
        if(!any(names(ellipsis)=="main")){
            ellipsis$main <- "Actuals vs Fitted";
        }
        # If type and ylab are not provided, set them...
        if(!any(names(ellipsis)=="type")){
            ellipsis$type <- "p";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- "Actuals";
        }
        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Fitted";
        }
        # xlim and ylim
        if(!any(names(ellipsis)=="xlim")){
            ellipsis$xlim <- range(c(ellipsis$x,ellipsis$y));
        }
        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- range(c(ellipsis$x,ellipsis$y));
        }

        # Start plotting
        do.call(plot,ellipsis);
        abline(a=0,b=1,col="grey",lwd=2,lty=2)
        if(lowess){
            lines(lowess(ellipsis$x, ellipsis$y), col="red");
        }
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

        if(is.occurrence(x$occurrence)){
            ellipsis$x <- ellipsis$x[actuals(x$occurrence)!=0];
            ellipsis$y <- ellipsis$y[actuals(x$occurrence)!=0];
        }

        # Remove NAs
        if(any(is.na(ellipsis$x))){
            ellipsis$x <- ellipsis$x[!is.na(ellipsis$x)];
            ellipsis$y <- ellipsis$y[!is.na(ellipsis$y)];
        }

        # Main, labs etc
        if(!any(names(ellipsis)=="main")){
            if(any(x$distribution==c("dinvgauss","dlnorm","dllaplace","dls"))){
                ellipsis$main <- paste0("log(",yName," Residuals) vs Fitted");
            }
            else{
                ellipsis$main <- paste0(yName," Residuals vs Fitted");
            }
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
                          "dlaplace"=,
                          "dllaplace"=qlaplace(c((1-level)/2, (1+level)/2), 0, 1),
                          "dalaplace"=qalaplace(c((1-level)/2, (1+level)/2), 0, 1, x$lambda),
                          "dlogis"=qlogis(c((1-level)/2, (1+level)/2), 0, 1),
                          "dt"=qt(c((1-level)/2, (1+level)/2), nobs(x)-nparam(x)),
                          "ds"=,
                          "dls"=qs(c((1-level)/2, (1+level)/2), 0, 1),
                          # In the next one, the scale is debiased, taking n-k into account
                          "dinvgauss"=qinvgauss(c((1-level)/2, (1+level)/2), mean=1,
                                                dispersion=x$scale * nobs(x) / (nobs(x)-nparam(x))),
                          "dlnorm"=qlnorm(c((1-level)/2, (1+level)/2), 0, 1),
                          qnorm(c((1-level)/2, (1+level)/2), 0, 1));
        # Analyse stuff in logarithms if the error is multiplicative
        if(any(x$distribution==c("dinvgauss","dlnorm"))){
            ellipsis$y[] <- log(ellipsis$y);
            zValues <- log(zValues);
        }
        else if(any(x$distribution==c("dllaplace","dls"))){
            ellipsis$y[] <- log(ellipsis$y);
        }
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
        abline(h=0, col="grey", lty=2);
        polygon(c(xRange,rev(xRange)),c(zValues[1],zValues[1],zValues[2],zValues[2]),
                col="lightgrey", border=NA, density=10);
        abline(h=zValues, col="red", lty=2);
        if(length(outliers)>0){
            points(ellipsis$x[outliers], ellipsis$y[outliers], pch=16);
            text(ellipsis$x[outliers], ellipsis$y[outliers], labels=outliers, pos=4);
        }
        if(lowess){
            lines(lowess(ellipsis$x[!is.na(ellipsis$y)], ellipsis$y[!is.na(ellipsis$y)]), col="red");
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
        ellipsis$y <- as.vector(residuals(x));
        if(any(x$distribution==c("dinvgauss","dlnorm","dllaplace","dls"))){
            ellipsis$y[] <- log(ellipsis$y);
        }
        if(type=="abs"){
            ellipsis$y[] <- abs(ellipsis$y);
        }
        else{
            ellipsis$y[] <- as.vector(ellipsis$y)^2;
        }

        if(is.occurrence(x$occurrence)){
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
                if(any(x$distribution==c("dinvgauss","dlnorm","dllaplace","dls"))){
                    ellipsis$main <- "|log(Residuals)| vs Fitted";
                }
                else{
                    ellipsis$main <- "|Residuals| vs Fitted";
                }
            }
            else{
                if(any(x$distribution==c("dinvgauss","dlnorm","dllaplace","dls"))){
                    ellipsis$main <- "log(Residuals)^2 vs Fitted";
                }
                else{
                    ellipsis$main <- "Residuals^2 vs Fitted";
                }
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
        abline(h=0, col="grey", lty=2);
        if(lowess){
            lines(lowess(ellipsis$x[!is.na(ellipsis$y)], ellipsis$y[!is.na(ellipsis$y)]), col="red");
        }
    }

    # 6. Q-Q with the specified distribution
    plot4 <- function(x, ...){
        ellipsis <- list(...);

        ellipsis$y <- residuals(x);
        if(is.occurrence(x$occurrence)){
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
                ellipsis$main <- "QQ plot of Normal distribution";
            }

            do.call(qqnorm, ellipsis);
            qqline(ellipsis$y);
        }
        else if(any(x$distribution=="dlnorm")){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ plot of Log Normal distribution";
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
        else if(x$distribution=="dllaplace"){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of Log Laplace distribution";
            }
            ellipsis$x <- exp(qlaplace(ppoints(500), mu=0, scale=x$scale));

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) exp(qlaplace(p, mu=0, scale=x$scale)));
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
        else if(x$distribution=="dls"){
            if(!any(names(ellipsis)=="main")){
                ellipsis$main <- "QQ-plot of Log S distribution";
            }
            ellipsis$x <- exp(qs(ppoints(500), mu=0, scale=x$scale));

            do.call(qqplot, ellipsis);
            qqline(ellipsis$y, distribution=function(p) exp(qs(p, mu=0, scale=x$scale)));
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

    # 7. Basic plot over time
    plot5 <- function(x, ...){
        yActuals <- actuals(x);
        if(!is.null(x$holdout)){
            yActuals <- ts(c(yActuals,x$holdout),start=start(yActuals),frequency=frequency(yActuals));
        }
        graphmaker(yActuals, x$forecast, fitted(x), main=x$model, legend=legend, parReset=FALSE, ...);
    }

    # 8 and 9. Standardised / Studentised residuals vs time
    plot6 <- function(x, type="rstandard", ...){

        ellipsis <- list(...);
        if(type=="rstandard"){
            ellipsis$x <- rstandard(x);
            yName <- "Standardised";
        }
        else{
            ellipsis$x <- rstudent(x);
            yName <- "Studentised";
        }

        if(is.occurrence(x$occurrence)){
            ellipsis$x <- ellipsis$x[actuals(x$occurrence)!=0];
        }

        if(!any(names(ellipsis)=="main")){
            ellipsis$main <- paste0(yName," Residuals vs Time");
        }

        if(!any(names(ellipsis)=="xlab")){
            ellipsis$xlab <- "Time";
        }
        if(!any(names(ellipsis)=="ylab")){
            ellipsis$ylab <- paste0(yName," Residuals");
        }

        # If type and ylab are not provided, set them...
        if(!any(names(ellipsis)=="type")){
            ellipsis$type <- "l";
        }

        zValues <- switch(x$distribution,
                          "dlaplace"=,
                          "dllaplace"=qlaplace(c((1-level)/2, (1+level)/2), 0, 1),
                          "dalaplace"=qalaplace(c((1-level)/2, (1+level)/2), 0, 1, x$lambda),
                          "dlogis"=qlogis(c((1-level)/2, (1+level)/2), 0, 1),
                          "dt"=qt(c((1-level)/2, (1+level)/2), nobs(x)-nparam(x)),
                          "ds"=,
                          "dls"=qs(c((1-level)/2, (1+level)/2), 0, 1),
                          # In the next one, the scale is debiased, taking n-k into account
                          "dinvgauss"=qinvgauss(c((1-level)/2, (1+level)/2), mean=1,
                                                dispersion=x$scale * nobs(x) / (nobs(x)-nparam(x))),
                          "dlnorm"=qlnorm(c((1-level)/2, (1+level)/2), 0, 1),
                          qnorm(c((1-level)/2, (1+level)/2), 0, 1));
        # Analyse stuff in logarithms if the error is multiplicative
        if(any(x$distribution==c("dinvgauss","dlnorm"))){
            ellipsis$x[] <- log(ellipsis$x);
            zValues <- log(zValues);
        }
        else if(any(x$distribution==c("dllaplace","dls"))){
            ellipsis$x[] <- log(ellipsis$x);
        }
        outliers <- which(ellipsis$x >zValues[2] | ellipsis$x <zValues[1]);


        if(!any(names(ellipsis)=="ylim")){
            ellipsis$ylim <- c(-max(abs(ellipsis$x)),max(abs(ellipsis$x)))*1.1;
        }

        if(legend){
            legendPosition <- "topright";
            ellipsis$ylim[2] <- ellipsis$ylim[2] + 0.2*diff(ellipsis$ylim);
            ellipsis$ylim[1] <- ellipsis$ylim[1] - 0.2*diff(ellipsis$ylim);
        }

        # Start plotting
        do.call(plot,ellipsis);
        if(length(outliers)>0){
            points(time(ellipsis$x)[outliers], ellipsis$x[outliers], pch=16);
            text(time(ellipsis$x)[outliers], ellipsis$x[outliers], labels=outliers, pos=4);
        }
        if(lowess){
            lines(lowess(c(1:length(ellipsis$x)),ellipsis$x), col="red");
        }
        abline(h=0, col="grey", lty=2);
        abline(h=zValues[1], col="red", lty=2);
        abline(h=zValues[2], col="red", lty=2);
        polygon(c(1:nobs(x), c(nobs(x):1)),
                c(rep(zValues[1],nobs(x)), rep(zValues[2],nobs(x))),
                col="lightgrey", border=NA, density=10);
        if(legend){
            legend(legendPosition,legend=c("Residuals",paste0(level*100,"% prediction interval")),
                   col=c("black","red"), lwd=rep(1,3), lty=c(1,1,2));
        }
    }

    # 10 and 11. ACF and PACF
    plot7 <- function(x, type="acf", ...){
        ellipsis <- list(...);

        if(!any(names(ellipsis)=="main")){
            if(type=="acf"){
                ellipsis$main <- "Autocorrelation Function of Residuals";
            }
            else{
                ellipsis$main <- "Partial Autocorrelation Function of Residuals";
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

    # 12. Plot of states
    plot8 <- function(x, ...){
        parDefault <- par(no.readonly = TRUE);
        if(any(unlist(gregexpr("C",x$model))==-1)){
            statesNames <- c("actuals",colnames(x$states),"residuals");
            x$states <- cbind(actuals(x),x$states,residuals(x));
            colnames(x$states) <- statesNames;
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
                if(ncol(x$states)<=5){
                    ellipsis$nc <- 1;
                }
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
        plot5(x, ...);
    }

    if(any(which==8)){
        plot6(x, ...);
    }

    if(any(which==9)){
        plot6(x, "rstudent", ...);
    }

    if(any(which==10)){
        plot7(x, type="acf", ...);
    }

    if(any(which==11)){
        plot7(x, type="pacf", ...);
    }

    if(any(which==12)){
        plot8(x, ...);
    }
}

#' @export
print.mes <- function(x, digits=4, ...){
    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"),2)," seconds"));
    cat(paste0("\nModel estimated: ",x$model));

    if(is.occurrence(x$occurrence)){
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
                             "general"="General",
                             "p"=,
                             "provided"="Provided by user");
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
                      "dllaplace" = "Log Laplace",
                      "dls" = "Log S",
                      # "dbcnorm" = paste0("Box-Cox Normal with lambda=",round(x$other$lambda,2)),
                      "dinvgauss" = "Inverse Gaussian"
    );
    if(is.occurrence(x$occurrence)){
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
       (any(x$loss==c("MSE","MSEh","MSCE","GPL")) & any(x$distribution==c("dnorm","dlnorm"))) ||
       (any(x$loss==c("aMSE","aMSEh","aMSCE","aGPL")) & any(x$distribution==c("dnorm","dlnorm"))) ||
       (any(x$loss==c("MAE","MAEh","MACE")) & any(x$distribution==c("dlaplace","dllaplace"))) ||
       (any(x$loss==c("HAM","HAMh","CHAM")) & any(x$distribution==c("ds","dls")))){
           ICs <- c(AIC(x),AICc(x),BIC(x),BICc(x));
           names(ICs) <- c("AIC","AICc","BIC","BICc");
           cat("\nInformation criteria:\n");
           print(round(ICs,digits));
    }
    else{
        cat("\nInformation criteria are unavailable for the chosen loss & distribution.\n");
    }

    # If there are accuracy measures, print them out
    if(!is.null(x$accuracy)){
        cat("\nForecast errors:\n");
        if(is.null(x$occurrence)){
            cat(paste(paste0("ME: ",round(x$accuracy["ME"],3)),
                      paste0("MAE: ",round(x$accuracy["MAE"],3)),
                      paste0("RMSE: ",round(sqrt(x$accuracy["MSE"]),3),"\n")
                      # paste0("Bias: ",round(x$accuracy["cbias"],3)*100,"%"),
                      ,sep="; "));
            cat(paste(paste0("sCE: ",round(x$accuracy["sCE"],5)*100,"%"),
                      paste0("sMAE: ",round(x$accuracy["sMAE"],5)*100,"%"),
                      paste0("sMSE: ",round(x$accuracy["sMSE"],5)*100,"%\n")
                ,sep="; "));
            cat(paste(paste0("MASE: ",round(x$accuracy["MASE"],3)),
                      paste0("RMSSE: ",round(x$accuracy["RMSSE"],3)),
                      paste0("rMAE: ",round(x$accuracy["rMAE"],3)),
                      paste0("rRMSE: ",round(x$accuracy["rRMSE"],3),"\n")
                      ,sep="; "));
        }
        else{
            cat(paste(paste0("Bias: ",round(x$accuracy["cbias"],5)*100,"%"),
                      paste0("sMSE: ",round(x$accuracy["sMSE"],5)*100,"%"),
                      paste0("rRMSE: ",round(x$accuracy["rRMSE"],3)),
                      paste0("sPIS: ",round(x$accuracy["sPIS"],5)*100,"%"),
                      paste0("sCE: ",round(x$accuracy["sCE"],5)*100,"%\n"),sep="; "));
        }
    }
}

#' @export
print.mesCombined <- function(x, digits=4, ...){
    cat(paste0("Time elapsed: ",round(as.numeric(x$timeElapsed,units="secs"),2)," seconds"));
    cat(paste0("\nModel estimated: ",x$model));

    cat(paste0("\n\nNumber of models combined: ", length(x$ICw)));
    cat(paste0("\nLoss function type: ",x$models[[1]]$loss));
    cat("\nSample size: "); cat(nobs(x));
    cat("\nAverage number of estimated parameters: "); cat(round(nparam(x),digits=digits));
    cat("\nAverage number of degrees of freedom: "); cat(round(nobs(x)-nparam(x),digits=digits));

    if(!is.null(x$accuracy)){
        cat("\n\nForecast errors:\n");
        if(is.null(x$occurrence)){
            cat(paste(paste0("ME: ",round(x$accuracy["ME"],3)),
                      paste0("MAE: ",round(x$accuracy["MAE"],3)),
                      paste0("RMSE: ",round(sqrt(x$accuracy["MSE"]),3),"\n")
                      # paste0("Bias: ",round(x$accuracy["cbias"],3)*100,"%"),
                      ,sep="; "));
            cat(paste(paste0("sCE: ",round(x$accuracy["sCE"],5)*100,"%"),
                      paste0("sMAE: ",round(x$accuracy["sMAE"],5)*100,"%"),
                      paste0("sMSE: ",round(x$accuracy["sMSE"],5)*100,"%\n")
                ,sep="; "));
            cat(paste(paste0("MASE: ",round(x$accuracy["MASE"],3)),
                      paste0("RMSSE: ",round(x$accuracy["RMSSE"],3)),
                      paste0("rMAE: ",round(x$accuracy["rMAE"],3)),
                      paste0("rRMSE: ",round(x$accuracy["rRMSE"],3),"\n")
                      ,sep="; "));
        }
        else{
            cat(paste(paste0("Bias: ",round(x$accuracy["cbias"],5)*100,"%"),
                      paste0("sMSE: ",round(x$accuracy["sMSE"],5)*100,"%"),
                      paste0("rRMSE: ",round(x$accuracy["rRMSE"],3)),
                      paste0("sPIS: ",round(x$accuracy["sPIS"],5)*100,"%"),
                      paste0("sCE: ",round(x$accuracy["sCE"],5)*100,"%\n"),sep="; "));
        }
    }
}

#### Coefficients ####
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

    return(mesReturn);
}

#' @export
coef.mes <- function(object, ...){
    return(object$B);
}


#' @importFrom stats sigma
#' @export
sigma.mes <- function(object, ...){
    df <- (nobs(object, all=FALSE)-nparam(object));
    # If the sample is too small, then use biased estimator
    if(df<=0){
        df[] <- nparam(object);
    }
    return(sqrt(switch(object$distribution,
                       "dnorm"=,
                       "dlogis"=,
                       "dlaplace"=,
                       "dt"=,
                       "ds"=,
                       "dalaplace"=sum(residuals(object)^2),
                       "dlnorm"=,
                       "dllaplace"=,
                       "dls"=sum(log(residuals(object))^2),
                       "dinvgauss"=sum((residuals(object)-1)^2))
                /df));
}

#' @export
summary.mes <- function(object, level=0.95, ...){
    ourReturn <- list(model=object$model,responseName=all.vars(formula(object))[1]);

    occurrence <- NULL;
    if(is.occurrence(object$occurrence)){
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
    parametersValues <- coef(object);
    if(is.null(parametersValues)){
        if(!is.null(object$xreg) && all(object$xregPersistence!=0)){
            parametersValues <- c(object$persistence,object$xregPersistence,object$initial,object$xregInitial);
        }
        else{
            parametersValues <- c(object$persistence,object$initial);
        }
        warning(paste0("Parameters are not available. You have probably provided them in the model, ",
                       "so there was nothing to estimate. We extracted smoothing parameters and initials."),
                call.=FALSE);
    }
    parametersConfint[,2:3] <- parametersValues + parametersConfint[,2:3];
    parametersTable <- cbind(parametersValues,parametersConfint);
    rownames(parametersTable) <- rownames(parametersConfint);
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
       (any(object$loss==c("MAE","MAEh","MACE")) & any(object$distribution==c("dlaplace","dllaplace"))) ||
       (any(object$loss==c("HAM","HAMh","CHAM")) & any(object$distribution==c("ds","dls")))){
        ICs <- c(AIC(object),AICc(object),BIC(object),BICc(object));
        names(ICs) <- c("AIC","AICc","BIC","BICc");
        ourReturn$ICs <- ICs;
    }
    return(structure(ourReturn, class="summary.mes"));
}

#' @export
summary.mesCombined <- function(object, ...){
    return(print.mesCombined(object, ...));
}

#' @export
print.summary.mes <- function(x, ...){
    ellipsis <- list(...);
    if(!any(names(ellipsis)=="digits")){
        digits <- 4;
    }
    else{
        digits <- ellipsis$digits;
    }

    cat(paste0("Model estimated: ",x$model));
    cat(paste0("\nResponse variable: ", paste0(x$responseName,collapse="")));

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
                      "dllaplace" = "Log Laplace",
                      "dls" = "Log S",
                      # "dbcnorm" = paste0("Box-Cox Normal with lambda=",round(x$other$lambda,2)),
                      "dinvgauss" = "Inverse Gaussian"
    );
    if(!is.null(x$occurrence)){
        distrib <- paste0("\nMixture of Bernoulli and ", distrib);
    }
    cat(paste0("\nDistribution used in the estimation: ", distrib));

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
       (any(x$loss==c("MAE","MAEh","MACE")) & any(x$distribution==c("dlaplace","dllaplace"))) ||
       (any(x$loss==c("HAM","HAMh","CHAM")) & any(x$distribution==c("ds","dls")))){
        cat("\nInformation criteria:\n");
        print(round(x$ICs,digits));
    }
    else{
        cat("\nInformation criteria are unavailable for the chosen loss & distribution.\n");
    }
}

#' @export
vcov.mes <- function(object, ...){
    # If the forecast is in numbers, then use its length as a horizon
    if(any(!is.na(object$forecast))){
        h <- length(object$forecast)
    }
    else{
        h <- 0;
    }
    modelReturn <- suppressWarnings(mes(actuals(object), h=h, model=object, FI=TRUE));
    vcovMatrix <- try(chol2inv(chol(modelReturn$FI)), silent=TRUE);
    if(inherits(vcovMatrix,"try-error")){
        vcovMatrix <- try(solve(modelReturn$FI, diag(ncol(modelReturn$FI)), tol=1e-20), silent=TRUE);
        if(inherits(vcovMatrix,"try-error")){
            warning(paste0("Sorry, but the hessian is singular, so we could not invert it.\n",
                           "We failed to produce the covariance matrix of parameters."),
                    call.=FALSE);
            vcovMatrix <- diag(1e+100,ncol(modelReturn$FI));
        }
    }
    colnames(vcovMatrix) <- rownames(vcovMatrix) <- colnames(modelReturn$FI);
    return(vcovMatrix);
}

#### Residuals and actuals functions ####

#' @importFrom greybox actuals
#' @export
actuals.mes <- function(object, all=TRUE, ...){
    if(all){
        return(object$y);
    }
    else{
        return(object$y[object$y!=0]);
    }
}

#' @export
nobs.mes <- function(object, ...){
    return(length(actuals(object, ...)));
}

#' @export
residuals.mes <- function(object, ...){
    return(switch(object$distribution,
                  "dlnorm"=,
                  "dllaplace"=,
                  "dls"=,
                  "dinvgauss"=switch(errorType(object),
                                     "A"=1+object$residuals/fitted(object),
                                     "M"=1+object$residuals),
                  "dnorm"=,
                  "dlogis"=,
                  "dlaplace"=,
                  "dt"=,
                  "ds"=,
                  "dalaplace"=,
                  object$residuals));
}

#' Multiple steps ahead forecast errors
#'
#' The function extracts 1 to h steps ahead forecast errors from the model.
#'
#' The errors correspond to the error term epsilon_t in the ETS models. Don't forget
#' that different models make different assumptions about epsilon_t and / or 1+epsilon_t.
#'
#' @template ssAuthor
#' @template ssKeywords
#'
#' @param object Model estimated using one of the forecasting functions.
#' @param h The forecasting horizon to use.
#' @param ... Currently nothing is accepted via ellipsis.
#' @return The matrix with observations in rows and h steps ahead values in columns.
#' So, the first row corresponds to the forecast produced from the 0th observation
#' from 1 to h steps ahead.
#' @seealso \link[stats]{residuals}, \link[stats]{rstandard}, \link[stats]{rstudent}
#' @examples
#'
#' x <- rnorm(100,0,1)
#' ourModel <- mes(x)
#' rmultistep(ourModel, h=13)
#'
#' @export rmultistep
rmultistep <- function(object, h=10, ...) UseMethod("rmultistep")

#' @export
rmultistep.default <- function(object, h=10, ...){
    return(NULL);
}

#' @export
rmultistep.mes <- function(object, h=10, ...){
    # Technical parameters
    lagsModelAll <- lags(object);
    componentsNumber <- length(object$persistence);
    componentsNumberSeasonal <- sum(lagsModelAll>1);
    lagsModelMax <- max(lagsModelAll);
    obsInSample <- nobs(object);

    # Model type
    model <- modelType(object);
    Etype <- errorType(object);
    Ttype <- substr(model,2,2);
    Stype <- substr(model,nchar(model),nchar(model));

    # Function returns the matrix with multi-step errors
    if(is.occurrence(object$occurrence)){
        ot <- matrix(actuals(object$occurrence),obsInSample,1);
    }
    else{
        ot <- matrix(1,obsInSample,1);
    }

    # Produce multi-step errors matrix
    return(ts(mesErrorerWrap(t(object$states), object$measurement, object$transition,
                             lagsModelAll, Etype, Ttype, Stype,
                             componentsNumber, componentsNumberSeasonal, h,
                             matrix(actuals(object),obsInSample,1), ot),
              start=start(actuals(object)), frequency=frequency(actuals(object))));
}

#' @importFrom stats rstandard
#' @export
rstandard.mes <- function(model, ...){
    obs <- nobs(model);
    df <- obs - nparam(model);
    errors <- residuals(model);
    # If this is an occurrence model, then only modify the non-zero obs
    # Also, if there are NAs in actuals, consider them as occurrence
    if(is.occurrence(model$occurrence)){
        residsToGo <- which(actuals(model$occurrence)!=0 & !is.na(actuals(model)));
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
    else if(model$distribution=="dls"){
        errors[] <- log(errors);
        return(exp((errors - mean(errors[residsToGo])) / (model$scale * obs / df)^2));
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
    # Also, if there are NAs in actuals, consider them as occurrence
    if(is.occurrence(model$occurrence)){
        residsToGo <- which(actuals(model$occurrence)!=0 & !is.na(actuals(model)));
    }
    else{
        residsToGo <- c(1:obs);
    }
    if(any(model$distribution==c("dt","dnorm"))){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / sqrt(sum(errors[-i]^2,na.rm=TRUE) / df);
        }
    }
    else if(model$distribution=="ds"){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sum(sqrt(abs(errors[-i])),na.rm=TRUE) / (2*df))^2;
        }
    }
    else if(model$distribution=="dls"){
        errors[] <- log(errors) - mean(log(errors));
        for(i in residsToGo){
            rstudentised[i] <- exp(errors[i] / (sum(sqrt(abs(errors[-i])),na.rm=TRUE) / (2*df))^2);
        }
    }
    else if(model$distribution=="dlaplace"){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sum(abs(errors[-i]),na.rm=TRUE) / df);
        }
    }
    else if(model$distribution=="dllaplace"){
        errors[] <- log(errors) - mean(log(errors));
        for(i in residsToGo){
            rstudentised[i] <- exp(errors[i] / (sum(abs(errors[-i]),na.rm=TRUE) / df));
        }
    }
    else if(model$distribution=="dalaplace"){
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sum(errors[-i] * (model$lambda - (errors[-i]<=0)*1),na.rm=TRUE) / df);
        }
    }
    else if(model$distribution=="dlogis"){
        errors[] <- errors - mean(errors);
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / (sqrt(sum(errors[-i]^2,na.rm=TRUE) / df) * sqrt(3) / pi);
        }
    }
    else if(model$distribution=="dlnorm"){
        errors[] <- log(errors) - mean(log(errors));
        for(i in residsToGo){
            rstudentised[i] <- exp(errors[i] / sqrt(sum(errors[-i]^2,na.rm=TRUE) / df));
        }
    }
    else if(model$distribution=="dinvgauss"){
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / mean(errors[residsToGo][-i],na.rm=TRUE);
        }
    }
    else{
        for(i in residsToGo){
            rstudentised[i] <- errors[i] / sqrt(sum(errors[-i]^2,na.rm=TRUE) / df);
        }
    }
    return(rstudentised);
}

#### Predict and forecast functions ####
#' @export
predict.mes <- function(object, newxreg=NULL, interval=c("none", "confidence", "prediction"),
                        level=0.95, side=c("both","upper","lower"), ...){

    interval <- match.arg(interval);
    obsInSample <- nobs(object);

    # Check if newxreg is provided
    if(!is.null(newxreg)){
        # If this is not a matrix / data.frame, then convert to one
        if(!is.data.frame(newxreg) && !is.matrix(newxreg)){
            newxreg <- as.data.frame(newxreg);
            colnames(newxreg) <- "xreg";
        }
        h <- nrow(newxreg);
        # If the newxreg is provided, then just do forecasts for that part
        if(any(interval==c("none","prediction","confidence"))){
            if(interval==c("prediction")){
                interval[] <- "simulated";
            }
            return(forecast(object, h=h, newxreg=newxreg,
                            interval=interval,
                            level=level, side=side, ...));
        }
    }
    else{
        # If there are no newxreg, then we need to produce fitted with / without interval
        if(interval=="none"){
            return(structure(list(mean=fitted(object), lower=NA, upper=NA, model=object,
                                  level=level, interval=interval, side=side),
                             class=c("mes.predict","mes.forecast")));
        }
        # Otherwise we do one-step-ahead prediction / confidence interval
        else{
            yForecast <- fitted(object);
        }
    }
    side <- match.arg(side);

    # Basic parameters
    model <- modelType(object);
    Etype <- errorType(object);

    # ts structure
    yForecastStart <- time(actuals(object))[obsInSample]+deltat(actuals(object));
    yFrequency <- frequency(actuals(object));

    # Extract variance and amend it in case of confidence interval
    s2 <- sigma(object)^2;
    if(interval=="confidence"){
        warning(paste0("Note that the ETS assumes that the initial level is known, ",
                       "so the confidence interval depends on smoothing parameters only."),
                call.=FALSE);
        s2 <- s2 * object$measurement[1:obsInSample,1:length(object$persistence),drop=FALSE] %*% object$persistence;
    }

    yUpper <- yLower <- yForecast;
    yUpper[] <- yLower[] <- NA;

    # If this is a mixture model, produce forecasts for the occurrence
    if(!is.null(object$occurrence)){
        occurrenceModel <- TRUE;
        pForecast <- fitted(object$occurrence);
    }
    else{
        occurrenceModel <- FALSE;
        pForecast <- rep(1, obsInSample);
    }

    # If this is an occurrence model, then take probability into account in the level.
    if(occurrenceModel && (interval=="prediction")){
        levelNew <- (level-(1-pForecast))/pForecast;
        levelNew[levelNew<0] <- 0;
    }
    else{
        levelNew <- level;
    }

    levelLow <- levelUp <- vector("numeric",obsInSample);
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

    #### Produce the intervals for the data ####
    if(object$distribution=="dnorm"){
        if(Etype=="A"){
            yLower[] <- qnorm(levelLow, 0, sqrt(s2));
            yUpper[] <- qnorm(levelUp, 0, sqrt(s2));
        }
        else{
            yLower[] <- qnorm(levelLow, 1, sqrt(s2));
            yUpper[] <- qnorm(levelUp, 1, sqrt(s2));
        }
    }
    else if(object$distribution=="dlogis"){
        if(Etype=="A"){
            yLower[] <- qlogis(levelLow, 0, sqrt(s2*3)/pi);
            yUpper[] <- qlogis(levelUp, 0, sqrt(s2*3)/pi);
        }
        else{
            yLower[] <- qlogis(levelLow, 1, sqrt(s2*3)/pi);
            yUpper[] <- qlogis(levelUp, 1, sqrt(s2*3)/pi);
        }
    }
    else if(object$distribution=="dlaplace"){
        if(Etype=="A"){
            yLower[] <- qlaplace(levelLow, 0, sqrt(s2/2));
            yUpper[] <- qlaplace(levelUp, 0, sqrt(s2/2));
        }
        else{
            yLower[] <- qlaplace(levelLow, 1, sqrt(s2/2));
            yUpper[] <- qlaplace(levelUp, 1, sqrt(s2/2));
        }
    }
    else if(object$distribution=="dt"){
        df <- nobs(object) - nparam(object);
        if(Etype=="A"){
            yLower[] <- sqrt(s2)*qt(levelLow, df);
            yUpper[] <- sqrt(s2)*qt(levelUp, df);
        }
        else{
            yLower[] <- (1 + sqrt(s2)*qt(levelLow, df));
            yUpper[] <- (1 + sqrt(s2)*qt(levelUp, df));
        }
    }
    else if(object$distribution=="ds"){
        if(Etype=="A"){
            yLower[] <- qs(levelLow, 0, (s2/120)^0.25);
            yUpper[] <- qs(levelUp, 0, (s2/120)^0.25);
        }
        else{
            yLower[] <- qs(levelLow, 1, (s2/120)^0.25);
            yUpper[] <- qs(levelUp, 1, (s2/120)^0.25);
        }
    }
    else if(object$distribution=="dalaplace"){
        lambda <- object$lambda;
        if(Etype=="A"){
            yLower[] <- qalaplace(levelLow, 0,
                                  sqrt(s2*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
            yUpper[] <- qalaplace(levelUp, 0,
                                  sqrt(s2*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
        }
        else{
            yLower[] <- qalaplace(levelLow, 1,
                                            sqrt(s2*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
            yUpper[] <- qalaplace(levelUp, 1,
                                            sqrt(s2*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
        }
    }
    else if(object$distribution=="dlnorm"){
        yLower[] <- qlnorm(levelLow, 0, sqrt(s2));
        yUpper[] <- qlnorm(levelUp, 0, sqrt(s2));
    }
    else if(object$distribution=="dllaplace"){
        yLower[] <- exp(qlaplace(levelLow, 0, sqrt(s2/2)));
        yUpper[] <- exp(qlaplace(levelUp, 0, sqrt(s2/2)));
    }
    else if(object$distribution=="dls"){
        yLower[] <- exp(qs(levelLow, 0, (s2/120)^0.25));
        yUpper[] <- exp(qs(levelUp, 0, (s2/120)^0.25));
    }
    else if(object$distribution=="dinvgauss"){
        yLower[] <- qinvgauss(levelLow, 1, dispersion=s2);
        yUpper[] <- qinvgauss(levelUp, 1, dispersion=s2);
    }

    #### Clean up the produced values for the interval ####
    # Make sensible values out of those weird quantiles
    if(Etype=="A"){
        yLower[levelLow==0] <- -Inf;
    }
    else{
        yLower[levelLow==0] <- 0;
    }
    yUpper[levelUp==1] <- Inf;

    # Substitute NAs and NaNs with zeroes
    if(any(is.nan(yLower)) || any(is.na(yLower))){
        yLower[is.nan(yLower)] <- 0;
        yLower[is.na(yLower)] <- 0;
    }
    if(any(is.nan(yUpper)) || any(is.na(yUpper))){
        yUpper[is.nan(yUpper)] <- 0;
        yUpper[is.na(yUpper)] <- 0;
    }

    if(Etype=="A"){
        yLower[] <- yForecast + yLower;
        yUpper[] <- yForecast + yUpper;
    }
    else{
        yLower[] <- yForecast * yLower;
        yUpper[] <- yForecast * yUpper;
    }

    return(structure(list(mean=yForecast, lower=yLower, upper=yUpper, model=object,
                          level=level, interval=interval, side=side),
                     class=c("mes.predict","mes.forecast")));
}

#' @export
plot.mes.predict <- function(x, ...){
    ellipsis <- list(...);
    if(is.null(ellipsis$ylim)){
        ellipsis$ylim <- range(c(actuals(x$model),x$mean,x$lower,x$upper),na.rm=TRUE);
    }
    ellipsis$x <- actuals(x$model);
    do.call(plot, ellipsis);
    lines(x$mean,col="purple",lwd=2,lty=2);
    if(x$interval!="none"){
        lines(x$lower,col="grey",lwd=3,lty=2);
        lines(x$upper,col="grey",lwd=3,lty=2);
    }
}

# Work in progress...
#' @param nsim Number of iterations to do in case of \code{interval="simulated"}.
#' @param occurrence The vector of occurrence variable (values in [0,1]).
#' @rdname forecast.smooth
#' @importFrom stats rnorm rlogis rt rlnorm qnorm qlogis qt qlnorm
#' @importFrom statmod rinvgauss qinvgauss
#' @importFrom greybox rlaplace rs ralaplace qlaplace qs qalaplace
#' @export
forecast.mes <- function(object, h=10, newxreg=NULL, occurrence=NULL,
                         interval=c("none", "simulated", "approximate", "semiparametric", "nonparametric", "confidence"),
                         level=0.95, side=c("both","upper","lower"), cumulative=FALSE, nsim=10000, ...){

    ellipsis <- list(...);

    interval <- match.arg(interval[1],c("none", "simulated", "approximate", "semiparametric",
                                        "nonparametric", "confidence", "parametric"));
    # If the horizon is zero, just construct fitted and potentially confidence interval thingy
    if(h<=0){
        if(all(interval!=c("none","confidence"))){
            interval[] <- "prediction";
        }
        return(predict(object, newxreg=newxreg,
                       interval=interval,
                       level=level, side=side, ...));
    }

    if(interval=="parametric"){
        warning("The parameter 'interval' does not accept 'parametric' anymore. We use 'approximate' value instead.",
                call.=FALSE);
        interval <- "approximate";
    }
    else if(interval=="confidence"){
        warning(paste0("Note that the ETS assumes that the initial level is known, ",
                       "so the confidence interval depends on smoothing parameters only."),
                call.=FALSE);
    }
    side <- match.arg(side);

    # Technical parameters
    lagsModelAll <- lags(object);
    componentsNumber <- length(lagsModelAll);
    componentsNumberSeasonal <- sum(lagsModelAll>1);
    lagsModelMax <- max(lagsModelAll);
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
    if(!is.null(object$xreg)){
        xregNumber <- ncol(object$xreg);
        if(is.null(newxreg)){
            warning("The newxreg is not provided. Predicting the explanatory variables based on what we have in-sample.",
                    call.=FALSE);
            newxreg <- matrix(NA,h,xregNumber);
            for(i in 1:xregNumber){
                newxreg[,i] <- mes(object$xreg[,i],h=h,silent=TRUE)$forecast;
            }
        }
        else{
            # If this is not a matrix / data.frame, then convert to one
            if(!is.data.frame(newxreg) && !is.matrix(newxreg)){
                newxreg <- as.data.frame(newxreg);
                colnames(newxreg) <- "xreg";
            }
            if(nrow(newxreg)<h){
                warning(paste0("The newxreg has ",nrow(newxreg)," observations, while ",h," are needed. ",
                               "Using the last available values as future ones."),
                        call.=FALSE);
                newnRows <- h-nrow(newxreg);
                xreg <- rbind(newxreg,matrix(rep(tail(newxreg,1),each=newnRows),newnRows,ncol(newxreg)));
            }
            else if(nrow(newxreg)>h){
                warning(paste0("The newxreg has ",nrow(newxreg)," observations, while only ",h," are needed. ",
                               "Using the last ",h," of them."),
                        call.=FALSE);
                xreg <- tail(newxreg,h);
            }
            else{
                xreg <- newxreg;
            }
            xregNames <- colnames(object$xreg);

            if(is.data.frame(xreg)){
                testFormula <- formula(object);
                testFormula[[2]] <- NULL;
                # Expand the variables and use only those that are in the model
                newxreg <- model.frame(testFormula, xreg);
                newxreg <- model.matrix(newxreg,data=newxreg)[,xregNames];
            }
            else{
                newxreg <- xreg[,xregNames];
            }
            rm(xreg);
        }

        componentsNumber[] <- componentsNumber - xregNumber;
        matWt[,componentsNumber+c(1:xregNumber)] <- newxreg;
        vecG <- matrix(c(object$persistence,object$xregPersistence), ncol=1);
    }
    else{
        vecG <- matrix(object$persistence, ncol=1);
        xregNumber <- 0;
    }
    matF <- object$transition;

    # Produce point forecasts
    mesForecast <- mesForecasterWrap(matVt, matWt, matF,
                                     lagsModelAll, Etype, Ttype, Stype,
                                     componentsNumber, componentsNumberSeasonal, h);

    #### Make safety checks
    # If there are NaN values
    if(any(is.nan(mesForecast))){
        mesForecast[is.nan(mesForecast)] <- 0;
    }
    # If there are negative values in the multiplicative model
    # if(any(c(Etype,Ttype,Stype)=="M") && any(mesForecast<=0)){
    #     mesForecast[mesForecast<=0] <- 0.01;
    # }

    # If this is a mixture model, produce forecasts for the occurrence
    if(is.occurrence(object$occurrence)){
        occurrenceModel <- TRUE;
        if(is.alm(object$occurrence)){
            pForecast <- forecast(object$occurrence,h=h,newdata=newxreg)$mean;
        }
        else{
            pForecast <- forecast(object$occurrence,h=h,newxreg=newxreg)$mean;
        }
    }
    else{
        occurrenceModel <- FALSE;
        # If this was provided occurrence, then use provided values
        if(!is.null(object$occurrence) && !is.null(object$occurrence$occurrence) &&
           (object$occurrence$occurrence=="provided")){
            if(!is.null(occurrence) && is.numeric(occurrence)){
                pForecast <- occurrence;
            }
            else{
                pForecast <- object$occurrence$forecast;
            }
            # Make sure that the values are of the correct length
            if(h<length(pForecast)){
                pForecast <- pForecast[1:h];
            }
            else if(h>length(pForecast)){
                pForecast <- c(pForecast,
                               rep(tail(pForecast,1),
                                   h-length(pForecast)));
            }
            else{
                pForecast <- pForecast;
            }
        }
        else{
            pForecast <- rep(1, h);
        }
    }

    # Cumulative forecasts have only one observation
    if(cumulative){
        yForecast <- yUpper <- yLower <- ts(vector("numeric", 1), start=yForecastStart, frequency=yFrequency);
        yForecast[] <- sum(mesForecast * pForecast);
    }
    else{
        yForecast <- yUpper <- yLower <- ts(vector("numeric", h), start=yForecastStart, frequency=yFrequency);
        yForecast[] <- mesForecast * pForecast;
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

    #### Simulated interval ####
    if(interval=="simulated"){
        arrVt <- array(NA, c(componentsNumber+xregNumber, h+lagsModelMax, nsim));
        arrVt[,1:lagsModelMax,] <- rep(matVt,nsim);
        sigmaValue <- sigma(object);
        matErrors <- matrix(switch(object$distribution,
                                   "dnorm"=rnorm(h*nsim, 0, sigmaValue),
                                   "dlogis"=rlogis(h*nsim, 0, sigmaValue*sqrt(3)/pi),
                                   "dlaplace"=rlaplace(h*nsim, 0, sigmaValue/2),
                                   "dt"=rt(h*nsim, obsInSample-nparam(object)),
                                   "ds"=rs(h*nsim, 0, (sigmaValue^2/120)^0.25),
                                   "dalaplace"=ralaplace(h*nsim, 0,
                                                         sqrt(sigmaValue^2*object$lambda^2*(1-object$lambda)^2/(object$lambda^2+(1-object$lambda)^2)),
                                                         object$lambda),
                                   "dlnorm"=rlnorm(h*nsim, 0, sigmaValue)-1,
                                   "dinvgauss"=rinvgauss(h*nsim, 1, dispersion=sigmaValue^2)-1,
                                   "dls"=exp(rs(h*nsim, 0, (sigmaValue^2/120)^0.25))-1,
                                   "dllaplace"=exp(rlaplace(h*nsim, 0, sigmaValue/2))-1
                                   ),
                            h,nsim);
        # This stuff is needed in order to produce adequate values for weird models
        EtypeModified <- Etype;
        if(Etype=="A" && any(object$distribution==c("dlnorm","dinvgauss","dls","dllaplace"))){
            EtypeModified[] <- "M";
        }

        # States, Errors, Ot, Transition, Measurement, Persistence
        ySimulated <- mesSimulatorwrap(arrVt, matErrors, matrix(rbinom(h*nsim, 1, pForecast), h, nsim),
                                       array(matF,c(dim(matF),nsim)), matWt,
                                       matrix(vecG, componentsNumber+xregNumber, nsim),
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
        # This step is needed in order to make intervals similar between the different methods
        if(Etype=="A"){
            yLower[] <- yLower - yForecast;
            yUpper[] <- yUpper - yForecast;
        }
        else{
            yLower[] <- yLower / yForecast;
            yUpper[] <- yUpper / yForecast;
        }
    }
    else{
        #### Approximate and confidence interval ####
        # Produce covatiance matrix and use it
        if(any(interval==c("approximate","confidence"))){
            s2 <- sigma(object)^2;
            # IG and Lnorm can use approximations from the multiplications
            if(any(object$distribution==c("dinvgauss","dlnorm","dls","dllaplace")) && Etype=="M"){
                vcovMulti <- mesVarAnal(lagsModelAll, h, matWt[1,,drop=FALSE], matF, vecG, s2);
                if(any(object$distribution==c("dlnorm","dls","dllaplace"))){
                    vcovMulti[] <- log(1+vcovMulti);
                }

                # The confidence interval relies on the assumption that initial level is known
                if(interval=="confidence"){
                    vcovMulti[] <- vcovMulti - s2;
                }

                # We don't do correct cumulatives in this case...
                if(cumulative){
                    vcovMulti <- sum(vcovMulti);
                }
            }
            else{
                vcovMulti <- covarAnal(lagsModelAll, h, matWt[1,,drop=FALSE], matF, vecG, s2);

                # The confidence interval relies on the assumption that initial level is known
                if(interval=="confidence"){
                    vcovMulti[] <- vcovMulti - s2;
                }

                # Do either the variance of sum, or a diagonal
                if(cumulative){
                    vcovMulti <- sum(vcovMulti);
                }
                else{
                    vcovMulti <- diag(vcovMulti);
                }
            }
        }
        #### Semiparametric and nonparametric interval ####
        # Extract multistep errors and calculate the covariance matrix
        else if(any(interval==c("semiparametric","nonparametric"))){
            if(h>1){
                mesErrors <- rmultistep(object, h=h);

                if(any(object$distribution==c("dinvgauss","dlnorm","dls","dllaplace")) && (Etype=="A")){
                    yFittedMatrix <- mesErrors;
                    for(i in 1:h){
                        yFittedMatrix[,i] <- fitted(object)[1:(obsInSample-h)+i];
                    }
                    mesErrors[] <- mesErrors/yFittedMatrix;
                }

                if(interval=="semiparametric"){
                    # Do either the variance of sum, or a diagonal
                    if(cumulative){
                        vcovMulti <- sum(t(mesErrors) %*% mesErrors / (obsInSample-h));
                    }
                    else{
                        vcovMulti <- diag(t(mesErrors) %*% mesErrors / (obsInSample-h));
                    }
                }
                # For nonparametric and cumulative...
                else{
                    if(cumulative){
                        mesErrors <- matrix(apply(mesErrors, 2, sum),obsInSample-h,1);
                    }
                }
            }
            else{
                vcovMulti <- sigma(object)^2;
                mesErrors <- resid(object);
            }
        }
        # Calculate interval for approximate and semiparametric
        if(any(interval==c("approximate","confidence","semiparametric"))){
            if(object$distribution=="dnorm"){
                if(Etype=="A"){
                    yLower[] <- qnorm(levelLow, 0, sqrt(vcovMulti));
                    yUpper[] <- qnorm(levelUp, 0, sqrt(vcovMulti));
                }
                else{
                    yLower[] <- qnorm(levelLow, 1, sqrt(vcovMulti));
                    yUpper[] <- qnorm(levelUp, 1, sqrt(vcovMulti));
                }
            }
            else if(object$distribution=="dlogis"){
                if(Etype=="A"){
                    yLower[] <- qlogis(levelLow, 0, sqrt(vcovMulti*3)/pi);
                    yUpper[] <- qlogis(levelUp, 0, sqrt(vcovMulti*3)/pi);
                }
                else{
                    yLower[] <- qlogis(levelLow, 1, sqrt(vcovMulti*3)/pi);
                    yUpper[] <- qlogis(levelUp, 1, sqrt(vcovMulti*3)/pi);
                }
            }
            else if(object$distribution=="dlaplace"){
                if(Etype=="A"){
                    yLower[] <- qlaplace(levelLow, 0, sqrt(vcovMulti/2));
                    yUpper[] <- qlaplace(levelUp, 0, sqrt(vcovMulti/2));
                }
                else{
                    yLower[] <- qlaplace(levelLow, 1, sqrt(vcovMulti/2));
                    yUpper[] <- qlaplace(levelUp, 1, sqrt(vcovMulti/2));
                }
            }
            else if(object$distribution=="dt"){
                df <- nobs(object) - nparam(object);
                if(Etype=="A"){
                    yLower[] <- sqrt(vcovMulti)*qt(levelLow, df);
                    yUpper[] <- sqrt(vcovMulti)*qt(levelUp, df);
                }
                else{
                    yLower[] <- (1 + sqrt(vcovMulti)*qt(levelLow, df));
                    yUpper[] <- (1 + sqrt(vcovMulti)*qt(levelUp, df));
                }
            }
            else if(object$distribution=="ds"){
                if(Etype=="A"){
                    yLower[] <- qs(levelLow, 0, (vcovMulti/120)^0.25);
                    yUpper[] <- qs(levelUp, 0, (vcovMulti/120)^0.25);
                }
                else{
                    yLower[] <- qs(levelLow, 1, (vcovMulti/120)^0.25);
                    yUpper[] <- qs(levelUp, 1, (vcovMulti/120)^0.25);
                }
            }
            else if(object$distribution=="dalaplace"){
                lambda <- object$lambda;
                if(Etype=="A"){
                    yLower[] <- qalaplace(levelLow, 0,
                                          sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
                    yUpper[] <- qalaplace(levelUp, 0,
                                          sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
                }
                else{
                    yLower[] <- qalaplace(levelLow, 1,
                                                    sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
                    yUpper[] <- qalaplace(levelUp, 1,
                                                    sqrt(vcovMulti*lambda^2*(1-lambda)^2/(lambda^2+(1-lambda)^2)), lambda);
                }
            }
            else if(object$distribution=="dlnorm"){
                yLower[] <- qlnorm(levelLow, 0, sqrt(vcovMulti));
                yUpper[] <- qlnorm(levelUp, 0, sqrt(vcovMulti));
            }
            else if(object$distribution=="dllaplace"){
                yLower[] <- exp(qlaplace(levelLow, 0, sqrt(vcovMulti/2)));
                yUpper[] <- exp(qlaplace(levelUp, 0, sqrt(vcovMulti/2)));
            }
            else if(object$distribution=="dls"){
                yLower[] <- exp(qs(levelLow, 0, (vcovMulti/120)^0.25));
                yUpper[] <- exp(qs(levelUp, 0, (vcovMulti/120)^0.25));
            }
            else if(object$distribution=="dinvgauss"){
                yLower[] <- qinvgauss(levelLow, 1, dispersion=vcovMulti);
                yUpper[] <- qinvgauss(levelUp, 1, dispersion=vcovMulti);
            }
        }
        # Use Taylor & Bunn approach for the nonparametric ones
        else if(interval=="nonparametric"){
            if(h>1){
                # This is needed in order to see if quant regression can be used
                if(all(levelLow==unique(levelLow))){
                    levelLow <- unique(levelLow);
                }
                if(all(levelUp==unique(levelUp))){
                    levelUp <- unique(levelUp);
                }

                # Do quantile regression for h>1 and scalars for the level
                if(length(levelLow)==1 && length(levelUp)==1){
                    # Quantile regression function
                    intervalQuantile <- function(A, alpha){
                        ee[] <- mesErrors - (A[1]*xe^A[2]);
                        return((1-alpha)*sum(abs(ee[ee<0]))+alpha*sum(abs(ee[ee>=0])));
                    }

                    ee <- mesErrors;
                    xe <- matrix(c(1:h),nrow=obsInSample-h,ncol=h,byrow=TRUE);

                    # lower quantiles
                    A <- nlminb(rep(1,2),intervalQuantile,alpha=levelLow)$par;
                    yLower[] <- A[1]*c(1:h)^A[2];

                    # upper quantiles
                    A[] <- nlminb(rep(1,2),intervalQuantile,alpha=levelUp)$par;
                    yUpper[] <- A[1]*c(1:h)^A[2];
                }
                else{
                    if(cumulative){
                        yLower[] <- quantile(mesErrors,levelLow,type=7);
                        yUpper[] <- quantile(mesErrors,levelUp,type=7);
                    }
                    else{
                        for(i in 1:h){
                            yLower[i] <- quantile(mesErrors[,i],levelLow[i],na.rm=T,type=7);
                            yUpper[i] <- quantile(mesErrors[,i],levelUp[i],na.rm=T,type=7);
                        }
                    }
                }
            }
            else{
                yLower[] <- quantile(mesErrors,levelLow,type=7);
                yUpper[] <- quantile(mesErrors,levelUp,type=7);
            }
            if(Etype=="M"){
                yLower[] <- 1+yLower;
                yUpper[] <- 1+yUpper;
            }
        }
        else{
            yUpper[] <- yLower[] <- NA;
        }
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

        # Substitute NAs and NaNs with zeroes
        if(any(is.nan(yLower)) || any(is.na(yLower))){
            yLower[is.nan(yLower)] <- switch(Etype,"A"=0,"M"=1);
            yLower[is.na(yLower)] <- switch(Etype,"A"=0,"M"=1);
        }
        if(any(is.nan(yUpper)) || any(is.na(yUpper))){
            yUpper[is.nan(yUpper)] <- switch(Etype,"A"=0,"M"=1);
            yUpper[is.na(yUpper)] <- switch(Etype,"A"=0,"M"=1);
        }

        # Do intervals around the forecasts...
        if(Etype=="A"){
            yLower[] <- yForecast + yLower;
            yUpper[] <- yForecast + yUpper;
        }
        else{
            yLower[] <- yForecast*yLower;
            yUpper[] <- yForecast*yUpper;
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

#' @export
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

#### Other methods ####
#' @export
multicov.mes <- function(object, type=c("analytical","empirical","simulated"), ...){
    type <- match.arg(type);

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

    if(type=="analytical"){
        covarMat <- covarAnal(lagsModelAll, h, matWt[1,,drop=FALSE], matF, vecG, s2);
    }
    else if(type=="empirical"){
        mesErrors <- rmultistep(object, h=h);
        covarMat <- t(mesErrors) %*% mesErrors / (nobs(object) - h);
    }

    return(covarMat);
}

#' @export
pointLik.mes <- function(object, ...){
    distribution <- object$distribution;
    yInSample <- actuals(object);
    obsInSample <- nobs(object);
    if(is.occurrence(object$occurrence)){
        otLogical <- yInSample!=0;
        yFitted <- fitted(object) / fitted(object$occurrence);
    }
    else{
        otLogical <- rep(TRUE, obsInSample);
        yFitted <- fitted(object);
    }
    scale <- object$scale;
    Etype <- errorType(object);

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
                                   "dllaplace"=dlaplace(q=log(yInSample[otLogical]), mu=log(yFitted[otLogical]),
                                                        scale=scale, log=TRUE),
                                   "dls"=ds(q=log(yInSample[otLogical]), mu=log(yFitted[otLogical]),
                                            scale=scale, log=TRUE),
                                   "dinvgauss"=dinvgauss(x=yInSample[otLogical], mean=yFitted[otLogical],
                                                         dispersion=scale/yFitted[otLogical], log=TRUE))

    # If this is a mixture model, take the respective probabilities into account (differential entropy)
    if(is.occurrence(object$occurrence)){
        likValues[!otLogical] <- -switch(distribution,
                                         "dnorm" =,
                                         "dlnorm" = (log(sqrt(2*pi)*scale)+0.5),
                                         "dlogis" = 2,
                                         "dlaplace" =,
                                         "dllaplace" =,
                                         "dalaplace" = (1 + log(2*scale)),
                                         "dt" = ((scale+1)/2 * (digamma((scale+1)/2)-digamma(scale/2)) +
                                                     log(sqrt(scale) * beta(scale/2,0.5))),
                                         "ds" =,
                                         "dls" = (2 + 2*log(2*scale)),
                                         "dinvgauss" = (0.5*(log(pi/2)+1+log(scale))));

        likValues[] <- likValues + pointLik(object$occurrence);
    }
    likValues <- ts(likValues, start=start(yFitted), frequency=frequency(yFitted));

    return(likValues);
}

##### Other methods to implement #####
# accuracy.mes <- function(object, holdout, ...){}
# simulate.mes <- function(object, nsim=1, seed=NULL, obs=NULL, ...){}
