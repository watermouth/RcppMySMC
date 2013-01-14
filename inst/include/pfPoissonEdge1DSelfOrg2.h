#ifndef PFPOISSONEDGE1DSELFORG2_H
#define PFPOISSONEDGE1DSELFORG2_H

#include "smctc.h"
#include <vector>
using namespace std;

namespace PoissonEdge1DSelfOrg2 {
  class State{
    public:
      vector<int> lineField;
      double alpha, beta;
      vector<double> gamma;
      double lambda;
  };
  double logLikelihood(long lTime, const State & x);
  smc::particle<State> fInitialise(smc::rng *pRng);
  void fMove(long lTime, smc::particle<State> & pFrom, smc::rng *pRng);
  double integrand_mean_x(const State& x, void *);
  double integrand_var_x(const State& x, void* vmx);
  double integrand_mean_lambda(const State& x, void*);
  double integrand_var_lambda(const State& x, void* mean);
}

#endif
