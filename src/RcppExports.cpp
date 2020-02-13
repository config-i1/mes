// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mesFitterWrap
RcppExport SEXP mesFitterWrap(SEXP matVt, SEXP matWt, SEXP matF, SEXP vecG, SEXP lagsModelAll, SEXP Etype, SEXP Ttype, SEXP Stype, SEXP componentsNumber, SEXP componentsNumberSeasonal, SEXP yInSample, SEXP ot, SEXP backcasting);
RcppExport SEXP _mes_mesFitterWrap(SEXP matVtSEXP, SEXP matWtSEXP, SEXP matFSEXP, SEXP vecGSEXP, SEXP lagsModelAllSEXP, SEXP EtypeSEXP, SEXP TtypeSEXP, SEXP StypeSEXP, SEXP componentsNumberSEXP, SEXP componentsNumberSeasonalSEXP, SEXP yInSampleSEXP, SEXP otSEXP, SEXP backcastingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matVt(matVtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type matWt(matWtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type matF(matFSEXP);
    Rcpp::traits::input_parameter< SEXP >::type vecG(vecGSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lagsModelAll(lagsModelAllSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Etype(EtypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Ttype(TtypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Stype(StypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type componentsNumber(componentsNumberSEXP);
    Rcpp::traits::input_parameter< SEXP >::type componentsNumberSeasonal(componentsNumberSeasonalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yInSample(yInSampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ot(otSEXP);
    Rcpp::traits::input_parameter< SEXP >::type backcasting(backcastingSEXP);
    rcpp_result_gen = Rcpp::wrap(mesFitterWrap(matVt, matWt, matF, vecG, lagsModelAll, Etype, Ttype, Stype, componentsNumber, componentsNumberSeasonal, yInSample, ot, backcasting));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mes_mesFitterWrap", (DL_FUNC) &_mes_mesFitterWrap, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_mes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
