#include "smctc.h"
#include "pfPoissonEdge1DSelfOrg2.h"
#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

namespace PoissonEdge1DSelfOrg2 {
  // ssm parameters
  double alpha0, beta0, gamma0, lambda0, lambdaMax;
  double sigma_alpha, sigma_gamma;
  ///The observations
  Rcpp::NumericVector y;
}

using namespace std;
using namespace Rcpp;
using namespace PoissonEdge1DSelfOrg2;

// [[Rcpp::export]]
SEXP pfPoissonEdge1DSelfOrg2(NumericVector observations, int numOfParticles, double sigma_alpha, double alpha0, double beta0,
                            double sigma_gamma, double gamma0,
                            double lambda0, double lambdaMax)
{
  BEGIN_RCPP
  // initial dummy observation
  y = observations;
  long lIterates = y.size();
  PoissonEdge1DSelfOrg2::sigma_alpha = sigma_alpha;
  PoissonEdge1DSelfOrg2::alpha0 = alpha0;
  PoissonEdge1DSelfOrg2::beta0  = beta0;
  PoissonEdge1DSelfOrg2::sigma_gamma = sigma_gamma;
  PoissonEdge1DSelfOrg2::gamma0 = gamma0;
  PoissonEdge1DSelfOrg2::lambda0 = lambda0;
  PoissonEdge1DSelfOrg2::lambdaMax= lambdaMax;
  //Initialise and run the sampler
  smc::sampler<State> Sampler(numOfParticles, SMC_HISTORY_NONE);  
  smc::moveset<State> Moveset(fInitialise, fMove, NULL);
  
  Sampler.SetResampleParams(SMC_RESAMPLE_MULTINOMIAL, 1.01 * numOfParticles);
  Sampler.SetMoveSet(Moveset);
  Sampler.Initialise();
  
  Rcpp::NumericVector meanEdge = Rcpp::NumericVector(lIterates);
  Rcpp::NumericVector meanLambda = Rcpp::NumericVector(lIterates);
  Rcpp::NumericVector varLambda = Rcpp::NumericVector(lIterates);
  vector<NumericVector> particleValues_lambda = vector<NumericVector>(lIterates);
  vector<NumericVector> particleValues_weight= vector<NumericVector>(lIterates);
  State tempState;
  for(int n=0 ; n < lIterates ; ++n) {
    if(n > 0) Sampler.Iterate();
    particleValues_lambda[n] = NumericVector(numOfParticles);
    particleValues_weight[n] = NumericVector(numOfParticles);
    for(int i=0; i<numOfParticles; i++){
      tempState = Sampler.GetParticleValue(i);
      (particleValues_lambda[n])[i] = tempState.lambda;
      (particleValues_weight[n])[i] = Sampler.GetParticleLogWeight(i);
    }
    meanEdge(n)   = Sampler.Integrate(integrand_mean_x,NULL);
    meanLambda(n) = Sampler.Integrate(integrand_mean_lambda,NULL);
    varLambda(n)  = Sampler.Integrate(integrand_var_lambda, (void*)&meanLambda(n));
  }
  return Rcpp::List::create(Rcpp::_["mean"] = meanEdge, Rcpp::_["meanLambda"] = meanLambda, Rcpp::_["varLambda"] = varLambda,
                            Rcpp::_["particles"] = particleValues_lambda,
                            Rcpp::_["pw"] = particleValues_weight);
  END_RCPP
}

namespace PoissonEdge1DSelfOrg2 {
  ///The function corresponding to the log likelihood at specified time and position (up to normalisation)
  
  ///  \param lTime The current time (i.e. the index of the current distribution)
  ///  \param X     The state to consider 
  double logLikelihood(long lTime, const State& x)
  {
    double ll = 0.0;
    double logYCum=0.0;
//    if(y[int(lTime)] > 0){
//      for(unsigned int i = 1; i < y[int(lTime)]; i++){
//        logYCum += log(i);
//      }
//    }
    if(x.lambda == 0){
      ll -= logYCum;
    } else {
      ll += y[int(lTime)] * log(x.lambda) - x.lambda - logYCum;
    }
//    ll += 0.5 * x.alpha * (1 - x.lineField[0]) * (x.gamma[0] * x.gamma[0]); // This is wrong idea...
    return ll;
  }
  
  ///A function to initialise particles
  
  /// \param pRng A pointer to the random number generator which is to be used
  smc::particle<State> fInitialise(smc::rng *pRng)
  {
    State x;
    x.lineField = vector<int>(1);
    x.lineField[0] = 0;
    x.alpha = alpha0;
    x.beta  = beta0;
    x.gamma = vector<double>(2);
    x.gamma[0] = pRng->Normal(gamma0,2*sigma_gamma);
    x.gamma[1] = pRng->Normal(gamma0,2*sigma_gamma);
    x.lambda = pRng->UniformDiscrete(lambda0, lambdaMax);
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
    double normalizationConstant, probOfNoEdge, termAlphaGamma;
    termAlphaGamma = exp(0.5 * to->alpha * (to->gamma[0] * to->gamma[0]));
    if(lineFieldAtPrevTime == 0){
      // no edge at previous time
      normalizationConstant = termAlphaGamma + 1;
    } else {
      normalizationConstant = termAlphaGamma + exp(-to->beta);
    }
    probOfNoEdge = termAlphaGamma / normalizationConstant;
    // transition
    if(pRng->UniformS() < probOfNoEdge){
      to->lineField[0] = 0;
      to->lambda = pRng->Normal(to->lambda, 1.0 / sqrt(to->alpha));
    } else {
      to->lineField[0] = 1;
      to->lambda = pRng->UniformDiscrete(1, lambdaMax); // model dependent min, max values
    }
//    std::cout << "lambda=" << to->lambda << std::endl;
    
    to->alpha     += pRng->Normal(0, sigma_alpha);
    to->beta      += pRng->Normal(0, beta0 / 100.0);
    to->gamma[0]  += pRng->Normal(0, sigma_gamma);
    double ll = logLikelihood(lTime, *to);
//    std::cout << "loglikelihood=" << ll << std::endl;
    pFrom.AddToLogWeight(ll);
  }  
  
  double integrand_mean_x(const State& x, void *)
  {
    return x.lineField[0];
  }
  double integrand_mean_lambda(const State& x, void *)
  {
    return x.lambda;
  }
  double integrand_var_lambda(const State& x, void* vmx)
  {
    double* dmx = (double*)vmx;
    double d = (x.lambda - (*dmx));
    return d*d;
  }

}

  
  
  