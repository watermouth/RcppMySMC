#include "smctc.h"
#include "pfEdge1DSelfOrg.h"
#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>

namespace Edge1DSelfOrg {
  // ssm parameters
  double alpha0, beta0, gamma0; 
  ///The observations
  Rcpp::NumericVector y;
}

using namespace std;
using namespace Rcpp;
using namespace Edge1DSelfOrg;

// [[Rcpp::export]]
SEXP pfEdge1DSelfOrg(NumericVector observations, int numOfParticles, double alpha0, double beta0, double gamma0)
{
  BEGIN_RCPP
  // initial dummy observation
  y = observations;
  long lIterates = y.size();
  Edge1DSelfOrg::alpha0 = alpha0;
  Edge1DSelfOrg::beta0  = beta0;
  Edge1DSelfOrg::gamma0 = gamma0;
  
  //Initialise and run the sampler
  smc::sampler<State> Sampler(numOfParticles, SMC_HISTORY_NONE);  
  smc::moveset<State> Moveset(fInitialise, fMove, NULL);
  
  Sampler.SetResampleParams(SMC_RESAMPLE_MULTINOMIAL, 1.01 * numOfParticles);
  Sampler.SetMoveSet(Moveset);
  Sampler.Initialise();
  
  Rcpp::NumericVector resMean = Rcpp::NumericVector(lIterates);
  for(int n=0 ; n < lIterates ; ++n) {
    if(n > 0) Sampler.Iterate();
    resMean(n) = Sampler.Integrate(integrand_mean_x,NULL);
  }
  return Rcpp::List::create(Rcpp::_["mean"] = resMean);
  END_RCPP
}

namespace Edge1DSelfOrg {
  ///The function corresponding to the log likelihood at specified time and position (up to normalisation)
  
  ///  \param lTime The current time (i.e. the index of the current distribution)
  ///  \param X     The state to consider 
  double logLikelihood(long lTime, const State& x)
  {
    double diffy;
    if(lTime == 0) {
      diffy = 0.0;
    } else {
      diffy = y[int(lTime)] - y[int(lTime)-1];
    }
    return - 0.5 * x.alpha * (1 - x.lineField[0]) * (diffy * diffy - x.gamma * x.gamma);
  }
  
  ///A function to initialise particles
  
  /// \param pRng A pointer to the random number generator which is to be used
  smc::particle<State> fInitialise(smc::rng *pRng)
  {
    State x;
    x.lineField = vector<int>(2);
    x.lineField[0] = 0;
    x.lineField[1] = 0;
    x.alpha = alpha0;
    x.beta  = beta0;
    x.gamma = gamma0;
    return smc::particle<State>(x,logLikelihood(0,x));
  }
  
  ///The proposal function.
  
  ///\param lTime The sampler iteration.
  ///\param pFrom The particle to move.
  ///\param pRng  A random number generator.
  void fMove(long lTime, smc::particle<State> & pFrom, smc::rng *pRng)
  {
    State *to = pFrom.GetValuePointer();
    int lineFieldAtPrevTime = to->lineField[0];
    // probability of a value which lineField0 takes is proportional to exp(-beta*KronekkerDelta(lineField0 + lineField1, 2))
    double normalizationConstant, probOfNoEdge;
    if(lineFieldAtPrevTime == 0){
      // no edge at previous time
      normalizationConstant = 1 + 1;
    } else {
      normalizationConstant = 1 + exp(-to->beta);
    }
    probOfNoEdge = 1 / normalizationConstant;
    // transition
    to->lineField[1] = to->lineField[0];
    if(pRng->UniformS() < probOfNoEdge){
      to->lineField[0] = 0;
    } else {
      to->lineField[0] = 1;
    }
    to->alpha += pRng->Normal(0, alpha0 / 100.0);
    to->beta += pRng->Normal(0, beta0 / 100.0);
    to->gamma += pRng->Normal(0, gamma0 / 100.0);
    pFrom.AddToLogWeight(logLikelihood(lTime, *to));
  }  
  
  double integrand_mean_x(const State& x, void *)
  {
    return x.lineField[0];
  }
}

  