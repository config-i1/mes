#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* # Function returns value of w() -- y-fitted -- used in the measurement equation */
inline double wvalue(arma::vec const &vecVt, arma::rowvec const &rowvecW,
                     char const &E, char const &T, char const &S,
                     unsigned int const &nNonSeasonal, unsigned int const &nSeasonal,
                     unsigned int const &nComponents){
    // vecVt is a vector here!
    double yfit = 0;
    arma::mat vecYfit;

    switch(S){
    // ZZN
    case 'N':
        switch(T){
        case 'N':
            vecYfit = rowvecW * vecVt;
        case 'A':
            vecYfit = rowvecW.cols(0,1) * vecVt.rows(0,1);
            break;
        case 'M':
            vecYfit = exp(rowvecW * log(vecVt));
            break;
        }
        break;
        // ZZA
    case 'A':
        switch(T){
        case 'N':
        case 'A':
            vecYfit = rowvecW * vecVt;
            break;
        case 'M':
            vecYfit = exp(rowvecW.cols(0,1) * log(vecVt.rows(0,1))) + rowvecW.cols(2,2+nSeasonal-1) * vecVt.rows(2,2+nSeasonal-1);
            break;
        }
        break;
        // ZZM
    case 'M':
        switch(T){
        case 'N':
        case 'M':
            vecYfit = exp(rowvecW * log(vecVt));
            break;
        case 'A':
            vecYfit = rowvecW.cols(0,1) * vecVt.rows(0,1) * exp(rowvecW.cols(2,2+nSeasonal-1) * log(vecVt.rows(2,2+nSeasonal-1)));
            break;
        }
        break;
    }

    // Explanatory variables
    if(nComponents > (nNonSeasonal+nSeasonal)){
        // If error is additive, add explanatory variables. Otherwise multiply by exp(ax)
        switch(E){
        case 'A':
            yfit = as_scalar(vecYfit + rowvecW.cols(nNonSeasonal+nSeasonal-1,nComponents-1) * vecVt.rows(nNonSeasonal+nSeasonal-1,nComponents-1));
            break;
        case 'M':
            yfit = as_scalar(vecYfit * exp(rowvecW.cols(nNonSeasonal+nSeasonal-1,nComponents-1) * vecVt.rows(nNonSeasonal+nSeasonal-1,nComponents-1)));
            break;
        }
    }

    return yfit;
}
//
// /* # Function returns value of r() -- additive or multiplicative error -- used in the error term of measurement equation.
//  This is mainly needed by sim.ets */
// inline double rvalue(arma::vec const &vecVt, arma::rowvec const &rowvecW, char const &E, char const &T, char const &S,
//                      arma::rowvec const &rowvecXt, arma::vec const &vecAt){
//
//     switch(E){
//     // MZZ
//     case 'M':
//         return wvalue(vecVt, rowvecW, E, T, S);
//         break;
//         // AZZ
//     case 'A':
//     default:
//         return 1.0;
//     }
// }

/* # Function returns value of f() -- new states without the update -- used in the transition equation */
inline arma::vec fvalue(arma::vec const &matrixVt, arma::mat const &matrixF, char const T, char const S,
                        unsigned int const &nComponents){
    arma::vec matrixVtnew = matrixVt;

    switch(T){
    case 'N':
    case 'A':
        matrixVtnew = matrixF * matrixVt;
        break;
    case 'M':
        matrixVtnew.rows(0,1) = exp(matrixF.submat(0,0,1,1) * log(matrixVt.rows(0,1)));
        if(nComponents>2){
            // This is needed in order not to face log(-x)
            matrixVtnew.rows(2,nComponents-1) = matrixVt.rows(2,nComponents-1);
        }
        break;
    }

    return matrixVtnew;
}

/* # Function returns value of g() -- the update of states -- used in components estimation for the persistence */
inline arma::vec gvalue(arma::vec const &matrixVt, arma::mat const &matrixF, arma::mat const &rowvecW,
                        char const &E, char const &T, char const &S,
                        unsigned int const &nNonSeasonal, unsigned int const &nSeasonal,
                        unsigned int const &nComponents){
    arma::vec g(matrixVt.n_rows, arma::fill::ones);

    // AZZ
    switch(E){
    case 'A':
        // ANZ
        switch(T){
        case 'N':
            switch(S){
            case 'M':
                g(0) = 1 / as_scalar(rowvecW.cols(1,nSeasonal-1) * matrixVt.rows(1,nSeasonal-1));
                g.rows(1,nSeasonal) = 1 / rowvecW(0) * matrixVt(0);
                // // Explanatory variables
                // if(nComponents > (nSeasonal+nNonSeasonal)){
                //     /* g.rows(1,nSeasonal) = 1 / (1/g(1) +
                //      as_scalar(rowvecW.cols(nSeasonal,nComponents-1) * matrixVt.rows(nSeasonal,nComponents))); */
                //     g.rows(nSeasonal+nNonSeasonal,nComponents-1) = 1/matrixVt.rows(nSeasonal+nNonSeasonal,nComponents-1);
                // }
                break;
            }
            break;
        // AAZ
        case 'A':
            switch(S){
            case 'M':
                g.rows(0,1) = g.rows(0,1) / as_scalar(rowvecW.cols(2,nSeasonal-1) * matrixVt.rows(2,nSeasonal-1));
                g.rows(2,2+nSeasonal-1) = 1 / as_scalar(rowvecW.cols(0,1) * matrixVt.rows(0,1));
                // // Explanatory variables
                // if(nComponents > (nSeasonal+nNonSeasonal)){
                //     /*g.rows(2,2+nSeasonal-1) = 1 / (1/g(1) +
                //       as_scalar(rowvecW.cols(nSeasonal,nComponents) * matrixVt.rows(nSeasonal,nComponents)));*/
                //     g.rows(nSeasonal+nNonSeasonal,nComponents-1) = 1/matrixVt.rows(nSeasonal+nNonSeasonal,nComponents-1);
                // }
                break;
            }
            break;
        // AMZ
        case 'M':
            switch(S){
            case 'N':
            case 'A':
                g(1) = g(1) / matrixVt(0);
                break;
            case 'M':
                g(0) = g(0) / as_scalar(rowvecW.cols(2,nSeasonal-1) * matrixVt.rows(2,nSeasonal-1));
                g(1) = g(1) / (matrixVt(0) * as_scalar(rowvecW.cols(2,nSeasonal-1) * matrixVt.rows(2,nSeasonal-1)));
                g.rows(2,2+nSeasonal-1) = g.rows(2,2+nSeasonal-1) / as_scalar(exp(rowvecW.cols(0,1) * log(matrixVt.rows(0,1))));
                break;
            }
            break;
        }
        break;
    // MZZ
    case 'M':
        // MNZ
        switch(T){
        case 'N':
            switch(S){
            case 'N':
                g = matrixVt;
                break;
            case 'A':
                g.rows(0,nSeasonal) = as_scalar(rowvecW.cols(0,nSeasonal) * matrixVt.rows(0,nSeasonal));
                // Explanatory variables
                if(nComponents > (nSeasonal+nNonSeasonal)){
                    g.rows(nSeasonal+nNonSeasonal,nComponents-1) = matrixVt.rows(nSeasonal+nNonSeasonal,nComponents-1);
                }
                break;
            case 'M':
                g = matrixVt;
                break;
            }
            break;
        // MAZ
        case 'A':
            switch(S){
            case 'N':
            case 'A':
                g.rows(0,nSeasonal) = as_scalar(rowvecW.cols(0,nSeasonal) * matrixVt.rows(0,nSeasonal));
                // Explanatory variables
                if(nComponents > (nSeasonal+nNonSeasonal)){
                    g.rows(nSeasonal+nNonSeasonal,nComponents-1) = matrixVt.rows(nSeasonal+nNonSeasonal,nComponents-1);
                }
                break;
            case 'M':
                g.rows(0,1).fill(as_scalar(rowvecW.cols(0,1) * matrixVt.rows(0,1)));
                g.rows(2,nComponents-1) = matrixVt.rows(2,nComponents-1);
                break;
            }
            break;
            // MMZ
        case 'M':
            switch(S){
            case 'N':
                g = exp(matrixF * log(matrixVt));
                break;
            case 'A':
                g.rows(0,2).fill(as_scalar(exp(rowvecW.cols(0,1) * log(matrixVt.rows(0,1))) + rowvecW.cols(2,nSeasonal-1) * matrixVt.rows(2,nSeasonal-1)));
                g(1) = g(0) / matrixVt(0);
                // Explanatory variables
                if(nComponents > (nSeasonal+nNonSeasonal)){
                    g.rows(nSeasonal+nNonSeasonal,nComponents-1) = matrixVt.rows(nSeasonal+nNonSeasonal,nComponents-1);
                }
                break;
            case 'M':
                g = exp(matrixF * log(matrixVt));
                break;
            }
            break;
        }
        break;
    }


    // Explanatory variables. Needed in order to update the parameters
    if(nComponents > (nSeasonal+nNonSeasonal)){
        arma::rowvec rowvecWtxreg(1/rowvecW.cols(nSeasonal+nNonSeasonal,nComponents-1).t());
        rowvecWtxreg.rows(find_nonfinite(rowvecWtxreg)).fill(0);
        switch(E){
        case 'A':
            g.rows(nSeasonal+nNonSeasonal,nComponents-1) = g.rows(nSeasonal+nNonSeasonal,nComponents-1)*rowvecWtxreg;
            break;
        case 'M':
            g.rows(nSeasonal+nNonSeasonal,nComponents-1) = g.rows(nSeasonal+nNonSeasonal,nComponents-1)*rowvecWtxreg;
            break;
        }
    }

    return g;
}

// /* # Function is needed for the renormalisation of seasonal components. It should be done seasonal-wise.*/
// inline arma::mat normaliser(arma::mat Vt, int &obsall, unsigned int &maxlag, char const &S, char const &T){
//
//     unsigned int nComponents = Vt.n_rows;
//     double meanseason = 0;
//
//     switch(S){
//     case 'A':
//         meanseason = mean(Vt.row(nComponents-1));
//         Vt.row(nComponents-1) = Vt.row(nComponents-1) - meanseason;
//         switch(T){
//         case 'N':
//         case 'A':
//             Vt.row(0) = Vt.row(0) + meanseason;
//             break;
//         case 'M':
//             Vt.row(0) = Vt.row(0) + meanseason / Vt.row(1);
//             break;
//         }
//         break;
//     case 'M':
//         meanseason = exp(mean(log(Vt.row(nComponents-1))));
//         Vt.row(nComponents-1) = Vt.row(nComponents-1) / meanseason;
//         switch(T){
//         case 'N':
//         case 'M':
//             Vt.row(0) = Vt.row(0) / meanseason;
//             break;
//         case 'A':
//             Vt.row(0) = Vt.row(0) * meanseason;
//             Vt.row(1) = Vt.row(1) * meanseason;
//             break;
//         }
//         break;
//     }
//
//     return(Vt);
// }
