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
#' \item \code{LASSO} - use LASSO to shrink the parameters of the model;
#' \item \code{MSE} (Mean Squared Error),
#' \item \code{MAE} (Mean Absolute Error),
#' \item \code{HAM} (Half Absolute Moment),
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
#' Starting values of parameters canbe passed via \code{x0}, while the upper and lower
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
#'
#' @return Object of class "mes" is returned. It contains the list of the
#' following values:
#'
#' @seealso \code{\link[forecast]{ets}, \link[smooth]{es}}
#'
#' @examples
#'
#' # Model selection using a specified pool of models
#' ourModel <- mes(rnorm(100,100,10),model=c("ANN","AAM","AMdA"))
#'
#' \dontrun{summary(ourModel)}
#' \dontrun{forecast(ourModel)}
#' \dontrun{plot(forecast(ourModel))}
#'
#' @importFrom forecast forecast
#' @importFrom greybox dlaplace dalaplace ds stepwise alm
#' @importFrom smooth modelType
#' @importFrom stats frequency dnorm dlogis dt dlnorm
#' @importFrom statmod dinvgauss
#' @importFrom nloptr nloptr
#' @importFrom numDeriv hessian
#' @export mes
mes <- function(y, model="ZZZ", lags=c(1,1,frequency(y)),
                distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                               "dlnorm","dinvgauss"),
                loss=c("likelihood","LASSO","MSE","MAE","HAM","MSEh","TMSE","GTMSE","MSCE"),
                h=0, holdout=FALSE,
                persistence=NULL, phi=NULL, initial=c("optimal","backcasting"),
                occurrence=c("none","auto","fixed","general","odds-ratio","inverse-odds-ratio","direct"),
                ic=c("AICc","AIC","BIC","BICc"), bounds=c("traditional","admissible","none"),
                xreg=NULL, xregDo=c("use","select"), xregInitial=NULL, xregPersistence=0,
                silent=TRUE, ...){
    # Copyright (C) 2019 - Inf  Ivan Svetunkov
    # Methods to implement:
    # cvar() - conditional variance,
    # predict(), forecast(),
    # plot() with options of what to plot,
    # resid() et al.,
    # vcov(), confint(),
    #
    # Parameters that were moved to forecast() and predict() functions:
    # h=10, holdout=FALSE, cumulative=FALSE,
    # interval=c("none","parametric","likelihood","semiparametric","nonparametric","confidence"), level=0.95,

    # Start measuring the time of calculations
    startTime <- Sys.time();

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
        #     persistence <- model$persistence;
        #     initial <- model$initial;
        #     initialSeason <- model$initialSeason;
        #     if(any(model$iprob!=1)){
        #         occurrence <- "a";
        #     }
        # }
        lags <- model$lags;
        date <- model$date;
        distribution <- model$distribution;
        loss <- model$loss;
        persistence <- model$persistence;
        phi <- model$phi;
        initial <- model$initial;
        occurrence <- model$occurrence;
        ic <- model$ic;
        bounds <- model$bounds;
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
    }
    else if(is.character(model)){
        # Everything is okay
    }
    else{
        warning("A model of an unknown class was provided. Switching to 'ZZZ'.",call.=FALSE);
        model <- "ZZZ";
    }

    #### Check the parameters of the function and create variables based on them ####
    parametersChecker(y, model, lags, persistence, phi, initial,
                      distribution, loss, h, holdout, occurrence, ic, bounds,
                      xreg, xregDo, xregInitial, xregPersistence,
                      silent, ParentEnvironment=environment(), ...);

    #### The function creates the necessary matrices based on the model and provided parameters ####
    # This is needed in order to initialise the estimation
    creator <- function(Etype, Ttype, Stype, lagsModel, lagsModelMax,
                        obsStates, obsInSample, componentsNumber, componentsNames,
                        yInSample, persistence, persistenceEstimate, phi,
                        initialValue, initialEstimate,
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
        if(initialEstimate){
            # For the seasonal models
            if(Stype!="N"){
                yDecomposition <- msdecompose(yInSample, lags[lags!=1], type=Stype);
                j <- 1;
                # level
                matVt[j,1:lagsModelMax] <- yDecomposition$initial[1];
                j <- j+1;
                if(Ttype!="N"){
                    matVt[j,1:lagsModelMax] <- switch(Ttype,
                                                      "A" = mean(diff(yDecomposition$trend),na.rm=TRUE),
                                                      "M" = exp(mean(diff(log(yDecomposition$trend)),na.rm=TRUE)));
                    j <- j+1;
                }
                for(i in j:componentsNumber){
                    matVt[i,(lagsModelMax-lagsModel[i])+1:lagsModel[i]] <- yDecomposition$seasonal[[i]];
                }
            }
            # Non-seasonal models
            else{
                matVt[1,1] <- mean(yInSample);
                if(Ttype!="N"){
                    matVt[2,1] <- switch(Ttype,
                                         "A" = mean(diff(yInSample),na.rm=TRUE),
                                         "M" = exp(mean(diff(log(yInSample)),na.rm=TRUE)));
                }
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
        if(xregInitialsProvided){
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
            vecG[j+1:xregNumber] <- B[j+1:xregNumber];
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
            matVt[componentsNumber+1:xregNumber,1:lagsModelMax] <- B[j+1:xregNumber];
        }

        return(list(matVt=matVt, matWt=matWt, matF=matF, vecG=vecG));
    }

    #### The function initialises the vector B for ETS ####
    initialiser <- function(Etype, Ttype, Stype,
                            componentsNumber, lagsModel, lagsModelMax, matVt,
                            persistenceEstimate, phiEstimate, initialType,
                            xregInitialsEstimate, xregPersistenceEstimate, xregNumber){
        # Persistence of states, persistence of xreg, phi, initials, initials for xreg
        B <- vector("numeric",persistenceEstimate*componentsNumber+xregPersistenceEstimate*xregNumber+phiEstimate+
                        (initialType=="optimal")*sum(lagsModel)+xregInitialsEstimate*xregNumber);
        j <- 1;
        # Fill in persistence
        if(persistenceEstimate){
            B[j:componentsNumber] <- switch(Etype,
                                            "A"=c(0.3,0.2,0.1),
                                            "M"=c(0.1,0.05,0.01))[j:componentsNumber];
            j <- j+componentsNumber;
        }

        # Persistence if xreg is provided
        if(xregPersistenceEstimate){
            B[j+1:xregNumber] <- rep(switch(Etype,
                                            "A"=0.1,
                                            "M"=0.01),xregNumber);
            j <- j+xregNumber;
        }

        # Damping parameter
        if(phiEstimate){
            B[j] <- 0.95;
            j <- j+1;
        }

        # Initials
        if(initialType=="optimal"){
            i <- 1;
            B[j] <- matVt[i,lagsModelMax];
            j <- j+1;
            i <- i+1;
            if(Ttype!="N"){
                B[j] <- matVt[i,lagsModelMax];
                j <- j+1;
                i <- i+1;
            }
            if(Stype!="N"){
                for(k in i:componentsNumber){
                    B[j+0:(lagsModel[k]-1)] <- matVt[k,(lagsModelMax-lagsModel[k])+1:lagsModel[k]];
                    j <- j+lagsModel[k];
                }
            }
        }

        # Initials of the xreg
        if(xregInitialsEstimate){
            B[j+1:xregNumber] <- matVt[componentsNumber+1:xregNumber,lagsModelMax];
        }

        return(B);
    }

    ##### Cost Function for ETS #####
    CF <- function(B,
                   Etype, Ttype, Stype, yInSample,
                   ot, otLogical, occurrenceModel,
                   componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                   matVt, matWt, matF, vecG, componentsNumberSeasonal,
                   persistenceEstimate, phiEstimate, initialType,
                   xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                   xregNumber,
                   bounds, loss, distribution, h, multisteps, other){
        ##### !!! distribution should not be default here !!! ####

        # Fill in the matrices
        mesElements <- filler(B,
                              Ttype, Stype, componentsNumber, lagsModel, lagsModelMax,
                              matVt, matWt, matF, vecG,
                              persistenceEstimate, phiEstimate, initialType,
                              xregInitialsEstimate, xregPersistenceEstimate, xregNumber);

        if(bounds=="traditional"){
            if(any(mesElements$vecG>1) || any(mesElements$vecG<0)){
                return(1E+300);
            }
        }
        else if(bounds=="admissible"){
            # We check the condition only for the last row of matWt
            eigenValues <- eigen(mesElements$matF - mesElements$vecG %*% mesElements$matWt[obsInSample,],
                                 only.values=TRUE)$values;
            if(any(eigenValues>1+1E-50)){
                return(eigenValues*1E+100);
            }
        }

        mesFitted <- mesFitterWrap(mesElements$matVt, mesElements$matWt, mesElements$matF, mesElements$vecG,
                                   lagsModelAll, Etype, Ttype, Stype, componentsNumber, componentsNumberSeasonal,
                                   yInSample, ot, initialType=="backcasting");

        if(!multisteps){
            if(loss=="likelihood"){
                # Scale for different functions
                scale <- switch(distribution,
                                "dnorm"=sqrt(sum(mesFitted$errors[otLogical]^2)/obsInSample),
                                "dlogis"=sqrt(sum(mesFitted$errors^2)/obsInSample * 3 / pi^2),
                                "dlaplace"=sum(abs(mesFitted$errors))/obsInSample,
                                "dt"=abs(other),
                                "ds"=sum(sqrt(abs(mesFitted$errors[otLogical]))) / (obsInSample*2),
                                "dalaplace"=sum(mesFitted$errors[otLogical]*(other-(mesFitted$errors[otLogical]<=0)*1))/obsInSample,
                                "dlnorm"=sqrt(sum(log(1+mesFitted$errors[otLogical])^2)/obsInSample),
                                "dinvgauss"=sum((mesFitted$errors[otLogical])^2/(1+mesFitted$errors[otLogical]))/obsInSample/mesFitted$yfit[otLogical]);

                # Calculate the likelihood
                CFValue <- -sum(switch(distribution,
                                       "dnorm"=dnorm(x=yInSample[otLogical], mean=mesFitted$yfit[otLogical],
                                                     sd=scale, log=TRUE),
                                       "dlogis"=dlogis(x=yInSample[otLogical], location=mesFitted$yfit[otLogical],
                                                       scale=scale, log=TRUE),
                                       "dlaplace"=dlaplace(q=yInSample[otLogical], mu=mesFitted$yfit[otLogical],
                                                           scale=scale, log=TRUE),
                                       "dt"=dt(mesFitted$errors, df=scale, log=TRUE),
                                       "ds"=ds(q=yInSample[otLogical],mu=mesFitted$yfit[otLogical],
                                               scale=scale, log=TRUE),
                                       "dalaplace"=dalaplace(q=yInSample[otLogical], mu=mesFitted$yfit[otLogical],
                                                             scale=scale, alpha=other, log=TRUE),
                                       "dlnorm"=dlnorm(x=yInSample[otLogical], meanlog=mesFitted$yfit[otLogical],
                                                       sdlog=scale, log=TRUE),
                                       "dinvgauss"=dinvgauss(x=yInSample[otLogical], mean=mesFitted$yfit[otLogical],
                                                             dispersion=scale, log=TRUE)));

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
            else if(loss=="LASSO"){
                # "B" needs to be normalised! Multiply by means of variables
                CFValue <- mean(mesFitted$errors^2) + other * sum(abs(B));
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
    logLikES <- function(B,
                         Etype, Ttype, Stype, yInSample,
                         ot, otLogical, occurrenceModel, pFitted,
                         componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                         matVt, matWt, matF, vecG, componentsNumberSeasonal,
                         persistenceEstimate, phiEstimate, initialType,
                         xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                         xregNumber,
                         bounds, loss, distribution, h, multisteps, other){
        if(!multisteps){
            if(loss=="LASSO"){
                return(0);
            }
            else{
                distributionNew <- switch(loss,
                                          "MSE"="dnorm",
                                          "MAE"="dlaplace",
                                          "HAM"="ds",
                                          distribution)
                logLikReturn <- CF(B,
                                   Etype, Ttype, Stype, yInSample, ot, otLogical, occurrenceModel,
                                   componentsNumber, lagsModel, lagsModelAll, lagsModelMax,
                                   matVt, matWt, matF, vecG, componentsNumberSeasonal,
                                   persistenceEstimate, phiEstimate, initialType,
                                   xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                                   xregNumber,
                                   bounds, "likelihood", distribution, h, multisteps, other);

                # If this is an occurrence model, add the probabilities
                if(occurrenceModel){
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
    estimator <- function(Etype, Ttype, Stype, lagsModel, lagsModelAll, lagsModelMax,
                          obsStates, obsInSample, componentsNumber, componentsNames,
                          componentsNumberSeasonal,
                          yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                          initialType, initialValue, initialEstimate,
                          xregProvided, xregInitialsProvided, xregInitialsEstimate,
                          xregPersistence, xregPersistenceEstimate,
                          xregModel, xregData, xregNumber, xregNames,
                          ot, otLogical, occurrenceModel, pFitted,
                          bounds, loss, distribution, h, multisteps, other){

        # Create the matrices for the specific ETS model
        esCreated <- creator(Etype, Ttype, Stype, lagsModel, lagsModelMax,
                             obsStates, obsInSample, componentsNumber, componentsNames,
                             yInSample, persistence, persistenceEstimate, phi,
                             initialValue, initialEstimate,
                             xregProvided, xregInitialsProvided, xregPersistence,
                             xregModel, xregData, xregNumber, xregNames);

        if(is.null(x0)){
            # Create the vector of initials for the optimisation
            B <- initialiser(Etype, Ttype, Stype,
                             componentsNumber, lagsModel, lagsModelMax, esCreated$matVt,
                             persistenceEstimate, phiEstimate, initialType,
                             xregInitialsEstimate, xregPersistenceEstimate, xregNumber);
        }
        else{
            B <- x0;
        }

        # Parameters are chosen to speed up the optimisation process and have decent accuracy
        res <- nloptr(B, CF, lb=lb, ub=ub,
                      opts=list(algorithm=algorithm, xtol_rel=xtol_rel, maxeval=maxeval, maxtime=maxtime, print_level=print_level),
                      Etype=Etype, Ttype=Ttype, Stype=Stype, yInSample=yInSample,
                      ot=ot, otLogical=otLogical, occurrenceModel=occurrenceModel,
                      componentsNumber=componentsNumber, lagsModel=lagsModel, lagsModelAll=lagsModelAll, lagsModelMax=lagsModelMax,
                      matVt=esCreated$matVt, matWt=esCreated$matWt, matF=esCreated$matF, vecG=esCreated$vecG,
                      componentsNumberSeasonal=componentsNumberSeasonal,
                      persistenceEstimate=persistenceEstimate, phiEstimate=phiEstimate, initialType=initialType,
                      xregProvided=xregProvided, xregInitialsEstimate=xregInitialsEstimate,
                      xregPersistenceEstimate=xregPersistenceEstimate, xregNumber=xregNumber,
                      bounds=bounds, loss=loss, distribution=distribution, h=h, multisteps=multisteps, other=other);

        # Prepare the values to return
        B[] <- res$solution;
        CFValue <- res$objective;
        # In case of likelihood, we typically have one more parameter to estimate - scale
        nParamEstimated <- length(B) + loss=="likelihood";
        logLikESValue <- logLikES(B,
                                  Etype, Ttype, Stype, yInSample,
                                  ot, otLogical, occurrenceModel, pFitted,
                                  componentsNumber, lagsModel, lagsModelMax,
                                  matVt, matWt, matF, vecG, componentsNumberSeasonal,
                                  persistenceEstimate, phiEstimate, initialType,
                                  xregProvided, xregInitialsEstimate, xregPersistenceEstimate,
                                  xregNumber,
                                  bounds, loss, distribution, h, multisteps, other);

        return(list(B=B, CFValue=CFValue, nParamEstimated=nParamEstimated, logLikESValue=logLikESValue));
    }

    ##### This function uses residuals in order to determine the needed xreg #####
    XregSelector <- function(listToReturn){
    }

    ##### Either estimate the model or create a pool #####
    if(modelDo=="estimate"){
        # Deal with occurrence model
        if(occurrenceModel && !occurrenceModelProvided){
            oesModel <- oes(yInSample, model=model, initial=initial, occurrence=occurrence, ic=ic, h=h,
                            holdout=FALSE, bounds=bounds, xreg=xreg, xregDo=xregDo);
            pFitted[] <- fitted(oesModel);
            parametersNumber[1,3] <- nparam(oesModel);
        }

        # Estimate the parameters of the demand sizes model
        esEstimator <- estimator(Etype, Ttype, Stype, lagsModel, lagsModelAll, lagsModelMax,
                                 obsStates, obsInSample, componentsNumber, componentsNames,
                                 componentsNumberSeasonal,
                                 yInSample, persistence, persistenceEstimate, phi, phiEstimate,
                                 initialType, initialValue, initialEstimate,
                                 xregProvided, xregInitialsProvided, xregInitialsEstimate,
                                 xregPersistence, xregPersistenceEstimate,
                                 xregModel, xregData, xregNumber, xregNames,
                                 ot, otLogical, occurrenceModel, pFitted,
                                 bounds, loss, distribution, h, multisteps, other);
        list2env(esEstimator);
        parametersNumber[1,4] <- nParamEstimated;

        # Create the matrices for the specific ETS model
        esCreated <- creator(Etype, Ttype, Stype, lagsModel, lagsModelMax,
                             obsStates, obsInSample, componentsNumber, componentsNames,
                             yInSample, persistence, persistenceEstimate, phi,
                             initialValue, initialEstimate,
                             xregProvided, xregInitialsProvided, xregPersistence,
                             xregModel, xregData, xregNumber, xregNames);
        list2env(esCreated);

        # Fill in the matrices
        mesElements <- filler(B,
                              Ttype, Stype, componentsNumber, lagsModel, lagsModelMax,
                              matVt, matWt, matF, vecG,
                              persistenceEstimate, phiEstimate, initialType,
                              xregInitialsEstimate, xregPersistenceEstimate, xregNumber);

        # Fit the model to the data
        mesFitted <- mesFitterWrap(mesElements$matVt, mesElements$matWt, mesElements$matF, mesElements$vecG,
                                   lagsModelAll, Etype, Ttype, Stype, componentsNumber, componentsNumberSeasonal,
                                   yInSample, ot, initialType=="backcasting");

        errors <- mesFitted$errors;
        yFitted <- mesFitted$yfit;
        matVt[] <- mesFitted$matVt;
    }

    return(structure(list(model=modelName, formula=mesFormula, timeElapsed=Sys.time()-startTime,
                          y=yInSample, holdout=yHoldout, fitted=yFitted, residuals=errors,
                          states=matVt, persistence=persistence, phi=phi, transition=matF,
                          measurement=matw, initialType=initialType, initial=initialValue,
                          nParam=parametersNumber, occurrence=occurrenceModel, xreg=xreg,
                          xregInitial=xregInitial, xregPersistence=xregPersistence,
                          ICs=ICs, logLik=logLikESValue, lossValue=CFValue, loss=loss,
                          distribution=distribution, B=B), class=c("mes","smooth")));
}
