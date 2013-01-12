#ifndef PFEDGE1D_H
#define PFEDGE1D_H

#include "smctc.h"
#include <vector>
using namespace std;

namespace Edge1D {
  class State{
    public:
      vector<int> lineField;
  };
  double logLikelihood(long lTime, const State & x);
  smc::particle<State> fInitialise(smc::rng *pRng);
  void fMove(long lTime, smc::particle<State> & pFrom, smc::rng *pRng);
  double integrand_mean_x(const State& x, void *);
  double integrand_var_x(const State& x, void* vmx);
}

#endif
