#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include "mesGeneral.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


/* # Function returns multiplicative or additive error for scalar */
double errorf(double const &yact, double &yfit, char const &E){
    switch(E){
        case 'A':
            return yact - yfit;
        break;
        case 'M':
            if((yact==0) & (yfit==0)){
                return 0;
            }
            else if((yact!=0) & (yfit==0)){
                return R_PosInf;
            }
            else{
                return (yact - yfit) / yfit;
            }
        break;
    }
}

/* # Function is needed to estimate the correct error for ETS when multisteps model selection with r(matvt) is sorted out. */
arma::mat errorvf(arma::mat yact, arma::mat yfit, char const &E){
    if(E=='A'){
        return yact - yfit;
    }
    else{
        yfit.elem(find(yfit==0)).fill(1e-100);
        return (yact - yfit) / yfit;
    }
}

// # Fitter for univariate models
List mesFitter(arma::mat &matrixVt, arma::mat const &matrixWt, arma::mat const &matrixF, arma::vec const &vectorG,
               arma::uvec &lags, char const &E, char const &T, char const &S,
               unsigned int const &nNonSeasonal, unsigned int const &nSeasonal,
               arma::vec const &vectorYt, arma::vec const &vectorOt, bool const &backcast){
    /* # matrixVt should have a length of obs + lagsModelMax.
    * # matrixWt is a matrix with nrows = obs
    * # vecG should be a vector
    * # lags is a vector of lags
    */

    int obs = vectorYt.n_rows;
    int obsall = matrixVt.n_cols;
    unsigned int nComponents = matrixVt.n_rows;
    arma::uvec lagsInternal = lags;
    unsigned int lagsModelMax = max(lagsInternal);
    int lagslength = lagsInternal.n_rows;

    lagsInternal = lagsInternal * nComponents;

    for(int i=0; i<lagslength; i=i+1){
        lagsInternal(i) = lagsInternal(i) + (lagslength - i - 1);
    }

    arma::uvec lagrows(lagslength, arma::fill::zeros);

    // Fitted values and the residuals
    arma::vec vecYfit(obs, arma::fill::zeros);
    arma::vec vecErrors(obs, arma::fill::zeros);

    // Loop for the backcasting
    if(backcast){
        unsigned int nIterations = 2;
    }
    else{
        unsigned int nIterations = 1;
    }
    // for (unsigned int j=1; j<nIterations; j=j+1) {

    // Loop for the model construction
    for (unsigned int i=lagsModelMax; i<obs+lagsModelMax; i=i+1) {
        lagrows = i * nComponents - lagsInternal + nComponents - 1;

        /* # Measurement equation and the error term */
        vecYfit(i-lagsModelMax) = wvalue(matrixVt(lagrows), matrixWt.row(i-lagsModelMax), E, T, S,
                nNonSeasonal, nSeasonal, nComponents);

        // This is a failsafe for cases of ridiculously high and ridiculously low values
        if(vecYfit(i-lagsModelMax) > 1e+100){
            vecYfit(i-lagsModelMax) = vecYfit(i-lagsModelMax-1);
        }

        // If this is zero (intermittent), then set error to zero
        if(vectorOt(i-lagsModelMax)==0){
            vecErrors(i-lagsModelMax) = 0;
        }
        else{
            vecErrors(i-lagsModelMax) = errorf(vectorYt(i-lagsModelMax), vecYfit(i-lagsModelMax), E);
        }

/* # Transition equation */
        matrixVt.col(i) = fvalue(matrixVt(lagrows), matrixF, T, S, nComponents) +
                          gvalue(matrixVt(lagrows), matrixF, matrixWt.row(i-lagsModelMax), E, T, S,
                                 nNonSeasonal, nSeasonal, nComponents) % vectorG * vecErrors(i-lagsModelMax);

/* Failsafe for cases when unreasonable value for state vector was produced */
        if(!matrixVt.col(i).is_finite()){
            matrixVt.col(i) = matrixVt(lagrows);
        }
        if((S=='M') & (matrixVt(matrixVt.n_rows-1,i) <= 0)){
            matrixVt(matrixVt.n_rows-1,i) = arma::as_scalar(matrixVt(lagrows.row(matrixVt.n_rows-1)));
        }
        if(T=='M'){
            if((matrixVt(0,i) <= 0) | (matrixVt(1,i) <= 0)){
                matrixVt(0,i) = arma::as_scalar(matrixVt(lagrows.row(0)));
                matrixVt(1,i) = arma::as_scalar(matrixVt(lagrows.row(1)));
            }
        }
        if(any(matrixVt.col(i)>1e+100)){
            matrixVt.col(i) = matrixVt(lagrows);
        }

/* Renormalise components if the seasonal model is chosen */
        // if(S!='N'){
        //     if(double(i+1) / double(lagsModelMax) == double((i+1) / lagsModelMax)){
        //         matrixVt.cols(i-lagsModelMax+1,i) = normaliser(matrixVt.cols(i-lagsModelMax+1,i), obsall, lagsModelMax, S, T);
        //     }
        // }

/* # Transition equation for xreg */
    }

    for (int i=obs+lagsModelMax; i<obsall; i=i+1) {
        lagrows = i * nComponents - lagsInternal + nComponents - 1;
        matrixVt.col(i) = fvalue(matrixVt(lagrows), matrixF, T, S, nComponents);

/* Failsafe for cases when unreasonable value for state vector was produced */
        if(!matrixVt.col(i).is_finite()){
            matrixVt.col(i) = matrixVt(lagrows);
        }
        if((S=='M') & (matrixVt(matrixVt.n_rows-1,i) <= 0)){
            matrixVt(matrixVt.n_rows-1,i) = arma::as_scalar(matrixVt(lagrows.row(matrixVt.n_rows-1)));
        }
        if(T=='M'){
            if((matrixVt(0,i) <= 0) | (matrixVt(1,i) <= 0)){
                matrixVt(0,i) = arma::as_scalar(matrixVt(lagrows.row(0)));
                matrixVt(1,i) = arma::as_scalar(matrixVt(lagrows.row(1)));
            }
        }
    }

    return List::create(Named("matvt") = matrixVt.t(), Named("yfit") = vecYfit,
                        Named("errors") = vecErrors);
}

/* # Wrapper for fitter */
// [[Rcpp::export]]
RcppExport SEXP mesFitterwrap(SEXP matVt, SEXP matWt, SEXP matF, SEXP vecG,
                              SEXP lagsModelAll, SEXP Etype, SEXP Ttype, SEXP Stype,
                              SEXP componentsNumber, SEXP componentsNumberSeasonal,
                              SEXP yInSample, SEXP ot, SEXP backcasting){

    NumericMatrix matvt_n(matVt);
    arma::mat matrixVt(matvt_n.begin(), matvt_n.nrow(), matvt_n.ncol());

    NumericMatrix matWt_n(matWt);
    arma::mat matrixWt(matWt_n.begin(), matWt_n.nrow(), matWt_n.ncol(), false);

    NumericMatrix matF_n(matF);
    arma::mat matrixF(matF_n.begin(), matF_n.nrow(), matF_n.ncol(), false);

    NumericMatrix vecg_n(vecG);
    arma::vec vectorG(vecg_n.begin(), vecg_n.nrow(), false);

    IntegerVector lagsModel_n(lagsModelAll);
    arma::uvec lags = as<arma::uvec>(lagsModel_n);

    char E = as<char>(Etype);
    char T = as<char>(Ttype);
    char S = as<char>(Stype);

    unsigned int nNonSeasonal = as<int>(componentsNumber);
    unsigned int nSeasonal = as<int>(componentsNumberSeasonal);

    NumericMatrix yt_n(yInSample);
    arma::vec vectorYt(yt_n.begin(), yt_n.nrow(), false);

    NumericVector ot_n(ot);
    arma::vec vectorOt(ot_n.begin(), ot_n.size(), false);

    bool backcast = as<char>(backcasting);

    return wrap(mesFitter(matrixVt, matrixWt, matrixF, vectorG,
                          lags, E, T, S,
                          nNonSeasonal, nSeasonal,
                          vectorYt, vectorOt, backcast));
}
