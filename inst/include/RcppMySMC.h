#ifndef RcppMySMC_h
#define RcppMySMC_h

#include "RcppMySMC_RcppExports.h"
#include <Rcpp.h>

namespace RcppMySMC {

    using namespace Rcpp;

    inline double onlineFiltering(smc::sampler<State>& Sampler, NumericVector observations) {
        static double(*p_onlineFiltering)(smc::sampler<State>&,NumericVector) = NULL;
        if (p_onlineFiltering == NULL) {
            validateExported("double(*onlineFiltering)(smc::sampler<State>&,NumericVector)");
            p_onlineFiltering = Rcpp::GetCppCallable("RcppMySMC", "RcppMySMC_RcppExports", "onlineFiltering");
        }
        return p_onlineFiltering(Sampler,observations);
    }
}

#endif // RcppMySMC_h
