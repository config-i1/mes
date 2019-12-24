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
#' \item \link[greybox]{dbcnorm} - Box-Cox normal distribution,
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
#' @param date The vector of dates for the corresponding values of \code{y}.
#' @param persistence Persistence vector \eqn{g}, containing smoothing
#' parameters. If \code{NULL}, then estimated.
#' @param phi Value of damping parameter. If \code{NULL} then it is estimated.
#' Only applicable for damped-trend models.
#' @param initial Can be either character or a vector of initial states. If it
#' is character, then it can be \code{"optimal"}, meaning that the initial
#' states are optimised, or \code{"backcasting"}, meaning that the initials are
#' produced using backcasting procedure (advised for data with high frequency).
#' @param loss The type of Loss Function used in optimization. \code{loss} can
#' be:
#' \itemize{
#' \item \code{likelihood} - the model is estimated via the maximisation of the
#' likelihood of the function specified in \code{distribution};
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
#' @param distribution what density function to assume for the error term. The full
#' name of the distribution should be provided, starting with the letter "d" -
#' "density". The names align with the names of distribution functions in R.
#' For example, see \link[stats]{dnorm}.
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
#' model, \code{usual} - restricting the values with (0, 1) or \code{none} - no
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
#' @param fast if \code{TRUE}, then the function won't check whether
#' the provided vectors are correct and will use them directly in the model
#' construction.
#' @param ...  Other non-documented parameters. For example \code{FI=TRUE} will
#' make the function also produce Fisher Information matrix, which then can be
#' used to calculated variances of smoothing parameters and initial states of
#' the model. This is used in the \link[stats]{vcov} method.
#' Parameters \code{A}, \code{ALower} and \code{AUpper} can be passed via
#' ellipsis as well. In this case they will be used for optimisation. \code{A}
#' sets the initial values before the optimisation, \code{ALower} and
#' \code{AUpper} define lower and upper bounds for the search inside of the
#' specified \code{bounds}. These values should have the length equal
#' to the number of parameters to estimate.
#' You can also pass parameters to the optimiser: 1. \code{maxeval} - maximum
#' number of evaluations to carry out (default is 100); 2. \code{xtol_rel} -
#' the precision of the optimiser (the default is 1E-6); 3. \code{algorithm} -
#' the algorithm to use in optimisation (\code{"NLOPT_LN_SBPLX"} by default).
#' 4. \code{print_level} - the level of output for the optimiser (0 by default).
#' You can read more about these parameters in the documentation of
#' \link[nloptr]{nloptr} function.
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
#' @importFrom greybox dlaplace dalaplace ds dbcnorm stepwise
#' @importFrom smooth modelType
#' @importFrom stats frequency dnorm dlogis dt dlnorm
#' @importFrom statmod dinvgauss
#' @importFrom nloptr nloptr
#' @importFrom numDeriv hessian
#' @export mes
mes <- function(y, model="ZZZ", lags=c(1,1,frequency(y)), date=NULL,
                 persistence=NULL, phi=NULL, initial=c("optimal","backcasting"),
                 loss=c("likelihood","MSE","MAE","HAM","MSEh","TMSE","GTMSE","MSCE"),
                 distribution=c("default","dnorm","dlogis","dlaplace","dt","ds","dalaplace",
                                "dlnorm","dbcnorm","dinvgauss"),
                 occurrence=c("none","auto","fixed","general","odds-ratio","inverse-odds-ratio","direct"),
                 ic=c("AICc","AIC","BIC","BICc"), bounds=c("usual","admissible","none"),
                 xreg=NULL, xregDo=c("use","select"), xregInitial=NULL, xregPersistence=0,
                 silent=TRUE, fast=FALSE, ...){
    # Copyright (C) 2019 - Inf  Ivan Svetunkov
    # Methods to implement:
    # cvar() - conditional variance, predict(), forecast(), plot() with options of what to plot, resid() et al., vcov(), confint(),
    #
    # Parameters that were moved to forecast() and predict() functions:
    # h=10, holdout=FALSE, cumulative=FALSE,
    # interval=c("none","parametric","likelihood","semiparametric","nonparametric","confidence"), level=0.95,

    # Start measuring the time of calculations
    startTime <- Sys.time();

# If a previous model provided as a model, write down the variables
    if(is.mes(model) || is.mes.sim(model)){
        if(is.omes(model$occurrence)){
            occurrence <- model$occurrence;
        }
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
        phi <- model$phi;
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
    # else if(forecast::is.ets(model)){
    #     # Extract smoothing parameters
    #     i <- 1;
    #     persistence <- coef(model)[i];
    #     if(model$components[2]!="N"){
    #         i <- i+1;
    #         persistence <- c(persistence,coef(model)[i]);
    #         if(model$components[3]!="N"){
    #             i <- i+1;
    #             persistence <- c(persistence,coef(model)[i]);
    #         }
    #     }
    #     else{
    #         if(model$components[3]!="N"){
    #             i <- i+1;
    #             persistence <- c(persistence,coef(model)[i]);
    #         }
    #     }
    #
    #     # Damping parameter
    #     if(model$components[4]=="TRUE"){
    #         i <- i+1;
    #         phi <- coef(model)[i];
    #     }
    #
    #     # Initials
    #     i <- i+1;
    #     initial <- coef(model)[i];
    #     if(model$components[2]!="N"){
    #         i <- i+1;
    #         initial <- c(initial,coef(model)[i]);
    #     }
    #
    #     # Initials of seasonal component
    #     if(model$components[3]!="N"){
    #         if(model$components[2]!="N"){
    #             initialSeason <- rev(model$states[1,-c(1:2)]);
    #         }
    #         else{
    #             initialSeason <- rev(model$states[1,-c(1)]);
    #         }
    #     }
    #     model <- modelType(model);
    # }
    else if(is.character(model)){
        # Everything is okay
    }
    else{
        warning("A model of an unknown class was provided. Switching to 'ZZZ'.",call.=FALSE);
        model <- "ZZZ";
    }

}
