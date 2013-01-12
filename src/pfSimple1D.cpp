#include "smctc.h"
#include "pfSimple1D.h"
#include <Rcpp.h>
#include <cstdlib>
#include <cmath>

namespace Simple1D {
double std_x0 = 10;
double std_x  = sqrt(10.0);
double var_y  = 1.0;

///The observations
Rcpp::NumericVector y;
}

using namespace std;
using namespace Rcpp;
using namespace Simple1D;

// [[Rcpp::export]]
SEXP pfSimple1D(NumericVector observations, int numOfParticles, double std_x, double std_y)
{
  y = observations;
  long lIterates = y.size();
  Simple1D::std_x = std_x;
  Simple1D::var_y = std_y * std_y;

  //Initialise and run the sampler
  smc::sampler<double> Sampler(numOfParticles, SMC_HISTORY_NONE);  
  smc::moveset<double> Moveset(fInitialise, fMove, NULL);

  Sampler.SetResampleParams(SMC_RESAMPLE_MULTINOMIAL, 1.01 * numOfParticles);
  Sampler.SetMoveSet(Moveset);
  Sampler.Initialise();

  Rcpp::NumericVector resMean = Rcpp::NumericVector(lIterates);
  Rcpp::NumericVector resSD   = Rcpp::NumericVector(lIterates);
  for(int n=0 ; n < lIterates ; ++n) {
    if(n > 0) Sampler.Iterate();
    resMean(n) = Sampler.Integrate(integrand_mean_x,NULL);
    resSD(n)  = sqrt(Sampler.Integrate(integrand_var_x, (void*)&resMean(n)));      
  }
  return Rcpp::List::create(Rcpp::_["mean"] = resMean,
			    Rcpp::_["sd"] = resSD);
}

namespace Simple1D {
///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
double logLikelihood(long lTime, const double & x)
{
    return -0.5 * pow(y[int(lTime)] - x, 2) / var_y;
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<double> fInitialise(smc::rng *pRng)
{
  double x;
  
  x = pRng->Normal(0,std_x0);

  return smc::particle<double>(x,logLikelihood(0,x));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<double> & pFrom, smc::rng *pRng)
{
  double *to = pFrom.GetValuePointer();

  *to += pRng->Normal(0.0, std_x);

  pFrom.AddToLogWeight(logLikelihood(lTime, *to));
}


double integrand_mean_x(const double& x, void *)
{
  return x;
}

double integrand_var_x(const double& x, void* vmx)
{
  double* dmx = (double*)vmx;
  double d = (x - (*dmx));
  return d*d;
}
}
