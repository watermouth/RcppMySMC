#ifndef PFSIMPLE1D_H
#define PFSIMPLE1D_H

#include "smctc.h"
using namespace std;

namespace Simple1D {
  double logLikelihood(long lTime, const double & x);
  smc::particle<double> fInitialise(smc::rng *pRng);
  void fMove(long lTime, smc::particle<double> & pFrom, smc::rng *pRng);
  double integrand_mean_x(const double& x, void *);
  double integrand_var_x(const double& x, void* vmx);
}

#endif

