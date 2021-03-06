// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_tiltmod_RCPPEXPORTS_H_GEN_
#define RCPP_tiltmod_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace tiltmod {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("tiltmod", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("tiltmod", "_tiltmod_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in tiltmod");
            }
        }
    }

    inline Rcpp::List MBM(const NumericVector& x, const NumericVector& w = NumericVector::create(), const NumericVector& a0 = NumericVector::create(), const NumericVector& a1 = NumericVector::create(), const double& precision = 1e-6, const int& MaxIter = 10000) {
        typedef SEXP(*Ptr_MBM)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_MBM p_MBM = NULL;
        if (p_MBM == NULL) {
            validateSignature("Rcpp::List(*MBM)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const double&,const int&)");
            p_MBM = (Ptr_MBM)R_GetCCallable("tiltmod", "_tiltmod_MBM");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_MBM(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(w)), Shield<SEXP>(Rcpp::wrap(a0)), Shield<SEXP>(Rcpp::wrap(a1)), Shield<SEXP>(Rcpp::wrap(precision)), Shield<SEXP>(Rcpp::wrap(MaxIter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List UBMM(const NumericVector& x, const NumericVector& w = NumericVector::create(), const NumericVector& a = NumericVector::create(), const double& precision = 1e-8, const int& MaxIter = 10000) {
        typedef SEXP(*Ptr_UBMM)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_UBMM p_UBMM = NULL;
        if (p_UBMM == NULL) {
            validateSignature("Rcpp::List(*UBMM)(const NumericVector&,const NumericVector&,const NumericVector&,const double&,const int&)");
            p_UBMM = (Ptr_UBMM)R_GetCCallable("tiltmod", "_tiltmod_UBMM");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_UBMM(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(w)), Shield<SEXP>(Rcpp::wrap(a)), Shield<SEXP>(Rcpp::wrap(precision)), Shield<SEXP>(Rcpp::wrap(MaxIter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

}

#endif // RCPP_tiltmod_RCPPEXPORTS_H_GEN_
