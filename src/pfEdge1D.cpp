#include "smctc.h"
#include "pfEdge1D.h"
#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>

namespace Edge1D {
  // ssm parameters
  double alpha, beta, gamma; 
  ///The observations
  Rcpp::NumericVector y;
}

using namespace std;
using namespace Rcpp;
using namespace Edge1D;

// [[Rcpp::export]]
SEXP pfEdge1D(NumericVector observations, int numOfParticles, double alpha, double beta, double gamma)
{
  BEGIN_RCPP
  // initial dummy observation
  y = observations;
  long lIterates = y.size();
  Edge1D::alpha = alpha;
  Edge1D::beta  = beta;
  Edge1D::gamma = gamma;
  
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

namespace Edge1D {
  ///The function corresponding to the log likelihood at specified time and position (up to normalisation)
  
  ///  \param lTime The current time (i.e. the index of the current distribution)
  ///  \param X     The state to consider 
  double logLikelihood(long lTime, const State& x)
  {
//    std::string number;
//    std::stringstream strstream;
//    strstream << lTime;
//    strstream >> number;
    double diffy;
    if(lTime == 0) {
      diffy = 0.0;
    } else {
      diffy = y[int(lTime)] - y[int(lTime)-1];
    }
//    if(lTime == 0){
//      std::cout << "lTime=" << lTime << ", y[-1]=" << y[-1] << ",diffy=" << diffy << std::endl;
//    }
    return - 0.5 * alpha * (1 - x.lineField[0]) * (diffy * diffy - gamma * gamma);
  }
  
  ///A function to initialise particles
  
  /// \param pRng A pointer to the random number generator which is to be used
  smc::particle<State> fInitialise(smc::rng *pRng)
  {
    State x;
    x.lineField = vector<int>(2);
    x.lineField[0] = 0;
    x.lineField[1] = 0;
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
      normalizationConstant = 1 + exp(-beta);
    }
    probOfNoEdge = 1 / normalizationConstant;
    // transition
    to->lineField[1] = to->lineField[0];
    if(pRng->UniformS() < probOfNoEdge){
      to->lineField[0] = 0;
    } else {
      to->lineField[0] = 1;
    }
    pFrom.AddToLogWeight(logLikelihood(lTime, *to));
  }  
  
  double integrand_mean_x(const State& x, void *)
  {
    return x.lineField[0];
  }
}

  