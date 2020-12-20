// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/myTridge.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// mvrnormArma
arma::mat mvrnormArma(const int n, arma::vec mu, const int p, const double rho);
static SEXP _myTridge_mvrnormArma_try(SEXP nSEXP, SEXP muSEXP, SEXP pSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, p, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _myTridge_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP pSEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_myTridge_mvrnormArma_try(nSEXP, muSEXP, pSEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rpois_rcpp
arma::vec rpois_rcpp(arma::vec& lambda);
static SEXP _myTridge_rpois_rcpp_try(SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rpois_rcpp(lambda));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _myTridge_rpois_rcpp(SEXP lambdaSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_myTridge_rpois_rcpp_try(lambdaSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cpprbinom
arma::vec cpprbinom(int n, double size, arma::vec prob);
static SEXP _myTridge_cpprbinom_try(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(cpprbinom(n, size, prob));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _myTridge_cpprbinom(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_myTridge_cpprbinom_try(nSEXP, sizeSEXP, probSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// genDataList
List genDataList(const int n, const arma::vec& mu, int p, double rho, arma::vec& beta, const double SNR, const std::string Test_case);
static SEXP _myTridge_genDataList_try(SEXP nSEXP, SEXP muSEXP, SEXP pSEXP, SEXP rhoSEXP, SEXP betaSEXP, SEXP SNRSEXP, SEXP Test_caseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type SNR(SNRSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    rcpp_result_gen = Rcpp::wrap(genDataList(n, mu, p, rho, beta, SNR, Test_case));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _myTridge_genDataList(SEXP nSEXP, SEXP muSEXP, SEXP pSEXP, SEXP rhoSEXP, SEXP betaSEXP, SEXP SNRSEXP, SEXP Test_caseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_myTridge_genDataList_try(nSEXP, muSEXP, pSEXP, rhoSEXP, betaSEXP, SNRSEXP, Test_caseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bFunction
arma::vec bFunction(arma::mat X, arma::vec theta, const std::string Test_case);
RcppExport SEXP _myTridge_bFunction(SEXP XSEXP, SEXP thetaSEXP, SEXP Test_caseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    rcpp_result_gen = Rcpp::wrap(bFunction(X, theta, Test_case));
    return rcpp_result_gen;
END_RCPP
}
// MeanFunction
arma::vec MeanFunction(arma::mat X, arma::vec theta, const std::string Test_case);
RcppExport SEXP _myTridge_MeanFunction(SEXP XSEXP, SEXP thetaSEXP, SEXP Test_caseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanFunction(X, theta, Test_case));
    return rcpp_result_gen;
END_RCPP
}
// MeanPrime
arma::vec MeanPrime(arma::mat X, arma::vec theta, const std::string Test_case);
RcppExport SEXP _myTridge_MeanPrime(SEXP XSEXP, SEXP thetaSEXP, SEXP Test_caseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanPrime(X, theta, Test_case));
    return rcpp_result_gen;
END_RCPP
}
// ObjectiveFunction
double ObjectiveFunction(arma::vec theta, const std::string Test_case, arma::vec y, arma::mat X, double trex_c);
RcppExport SEXP _myTridge_ObjectiveFunction(SEXP thetaSEXP, SEXP Test_caseSEXP, SEXP ySEXP, SEXP XSEXP, SEXP trex_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type trex_c(trex_cSEXP);
    rcpp_result_gen = Rcpp::wrap(ObjectiveFunction(theta, Test_case, y, X, trex_c));
    return rcpp_result_gen;
END_RCPP
}
// Gradient
arma::vec Gradient(arma::vec theta, const std::string Test_case, arma::vec y, arma::mat X, double trex_c);
RcppExport SEXP _myTridge_Gradient(SEXP thetaSEXP, SEXP Test_caseSEXP, SEXP ySEXP, SEXP XSEXP, SEXP trex_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type trex_c(trex_cSEXP);
    rcpp_result_gen = Rcpp::wrap(Gradient(theta, Test_case, y, X, trex_c));
    return rcpp_result_gen;
END_RCPP
}
// GradientLs
arma::vec GradientLs(arma::vec theta, arma::mat X, arma::vec y, const std::string Test_case);
RcppExport SEXP _myTridge_GradientLs(SEXP thetaSEXP, SEXP XSEXP, SEXP ySEXP, SEXP Test_caseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    rcpp_result_gen = Rcpp::wrap(GradientLs(theta, X, y, Test_case));
    return rcpp_result_gen;
END_RCPP
}
// optim_ObLs
arma::vec optim_ObLs(arma::vec theta, arma::mat X, arma::vec y, const std::string Test_case);
RcppExport SEXP _myTridge_optim_ObLs(SEXP thetaSEXP, SEXP XSEXP, SEXP ySEXP, SEXP Test_caseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_ObLs(theta, X, y, Test_case));
    return rcpp_result_gen;
END_RCPP
}
// optim_Ridge
arma::vec optim_Ridge(arma::vec theta, arma::mat X, arma::vec y, const std::string Test_case, double r);
RcppExport SEXP _myTridge_optim_Ridge(SEXP thetaSEXP, SEXP XSEXP, SEXP ySEXP, SEXP Test_caseSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_Ridge(theta, X, y, Test_case, r));
    return rcpp_result_gen;
END_RCPP
}
// optim_ObFn
arma::vec optim_ObFn(arma::vec theta, arma::mat X, arma::vec y, const std::string Test_case, double trex_c);
RcppExport SEXP _myTridge_optim_ObFn(SEXP thetaSEXP, SEXP XSEXP, SEXP ySEXP, SEXP Test_caseSEXP, SEXP trex_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::string >::type Test_case(Test_caseSEXP);
    Rcpp::traits::input_parameter< double >::type trex_c(trex_cSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_ObFn(theta, X, y, Test_case, trex_c));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _myTridge_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*mvrnormArma)(const int,arma::vec,const int,const double)");
        signatures.insert("arma::vec(*rpois_rcpp)(arma::vec&)");
        signatures.insert("arma::vec(*cpprbinom)(int,double,arma::vec)");
        signatures.insert("List(*genDataList)(const int,const arma::vec&,int,double,arma::vec&,const double,const std::string)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _myTridge_RcppExport_registerCCallable() { 
    R_RegisterCCallable("myTridge", "_myTridge_mvrnormArma", (DL_FUNC)_myTridge_mvrnormArma_try);
    R_RegisterCCallable("myTridge", "_myTridge_rpois_rcpp", (DL_FUNC)_myTridge_rpois_rcpp_try);
    R_RegisterCCallable("myTridge", "_myTridge_cpprbinom", (DL_FUNC)_myTridge_cpprbinom_try);
    R_RegisterCCallable("myTridge", "_myTridge_genDataList", (DL_FUNC)_myTridge_genDataList_try);
    R_RegisterCCallable("myTridge", "_myTridge_RcppExport_validate", (DL_FUNC)_myTridge_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_myTridge_mvrnormArma", (DL_FUNC) &_myTridge_mvrnormArma, 4},
    {"_myTridge_rpois_rcpp", (DL_FUNC) &_myTridge_rpois_rcpp, 1},
    {"_myTridge_cpprbinom", (DL_FUNC) &_myTridge_cpprbinom, 3},
    {"_myTridge_genDataList", (DL_FUNC) &_myTridge_genDataList, 7},
    {"_myTridge_bFunction", (DL_FUNC) &_myTridge_bFunction, 3},
    {"_myTridge_MeanFunction", (DL_FUNC) &_myTridge_MeanFunction, 3},
    {"_myTridge_MeanPrime", (DL_FUNC) &_myTridge_MeanPrime, 3},
    {"_myTridge_ObjectiveFunction", (DL_FUNC) &_myTridge_ObjectiveFunction, 5},
    {"_myTridge_Gradient", (DL_FUNC) &_myTridge_Gradient, 5},
    {"_myTridge_GradientLs", (DL_FUNC) &_myTridge_GradientLs, 4},
    {"_myTridge_optim_ObLs", (DL_FUNC) &_myTridge_optim_ObLs, 4},
    {"_myTridge_optim_Ridge", (DL_FUNC) &_myTridge_optim_Ridge, 5},
    {"_myTridge_optim_ObFn", (DL_FUNC) &_myTridge_optim_ObFn, 5},
    {"_myTridge_RcppExport_registerCCallable", (DL_FUNC) &_myTridge_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_myTridge(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
