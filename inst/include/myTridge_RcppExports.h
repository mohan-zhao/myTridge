// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_myTridge_RCPPEXPORTS_H_GEN_
#define RCPP_myTridge_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace myTridge {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("myTridge", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("myTridge", "_myTridge_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in myTridge");
            }
        }
    }

    inline arma::mat mvrnormArma(const int n, arma::vec mu, const int p, const double rho) {
        typedef SEXP(*Ptr_mvrnormArma)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_mvrnormArma p_mvrnormArma = NULL;
        if (p_mvrnormArma == NULL) {
            validateSignature("arma::mat(*mvrnormArma)(const int,arma::vec,const int,const double)");
            p_mvrnormArma = (Ptr_mvrnormArma)R_GetCCallable("myTridge", "_myTridge_mvrnormArma");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mvrnormArma(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(mu)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(rho)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::vec rpois_rcpp(arma::vec& lambda) {
        typedef SEXP(*Ptr_rpois_rcpp)(SEXP);
        static Ptr_rpois_rcpp p_rpois_rcpp = NULL;
        if (p_rpois_rcpp == NULL) {
            validateSignature("arma::vec(*rpois_rcpp)(arma::vec&)");
            p_rpois_rcpp = (Ptr_rpois_rcpp)R_GetCCallable("myTridge", "_myTridge_rpois_rcpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rpois_rcpp(Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec cpprbinom(int n, double size, arma::vec prob) {
        typedef SEXP(*Ptr_cpprbinom)(SEXP,SEXP,SEXP);
        static Ptr_cpprbinom p_cpprbinom = NULL;
        if (p_cpprbinom == NULL) {
            validateSignature("arma::vec(*cpprbinom)(int,double,arma::vec)");
            p_cpprbinom = (Ptr_cpprbinom)R_GetCCallable("myTridge", "_myTridge_cpprbinom");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpprbinom(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(size)), Shield<SEXP>(Rcpp::wrap(prob)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline List genDataList(const int n, const arma::vec& mu, int p, double rho, arma::vec& beta, const double SNR, const std::string Test_case) {
        typedef SEXP(*Ptr_genDataList)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_genDataList p_genDataList = NULL;
        if (p_genDataList == NULL) {
            validateSignature("List(*genDataList)(const int,const arma::vec&,int,double,arma::vec&,const double,const std::string)");
            p_genDataList = (Ptr_genDataList)R_GetCCallable("myTridge", "_myTridge_genDataList");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_genDataList(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(mu)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(SNR)), Shield<SEXP>(Rcpp::wrap(Test_case)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_myTridge_RCPPEXPORTS_H_GEN_