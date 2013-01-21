#include "smctc.h"
#include "sampler.h"
#include "pfPoissonEdge1DSelfOrg.h"

using namespace std;
using namespace Rcpp;
using namespace PoissonEdge1DSelfOrg;

namespace PoissonEdge1DSelfOrg{
  extern NumericVector y;
  extern double alpha0, beta0, gamma0, lambda0, lambdaMax;
  extern double sigma_alpha;
  smc::sampler<State>* p_sampler;
}

double onlineFiltering(smc::sampler<State>& Sampler, double observation, size_t i_obs){
  y[i_obs] = observation;
  Sampler.Iterate();
  return(Sampler.Integrate(integrand_mean_lambda,NULL));
}

void initialize(double initialObservation, int numOfParticles, double sigma_alpha, double alpha0, double beta0,
                            double lambda0, double lambdaMax)
{
  PoissonEdge1DSelfOrg::sigma_alpha = sigma_alpha;
  PoissonEdge1DSelfOrg::alpha0 = alpha0;
  PoissonEdge1DSelfOrg::beta0  = beta0;
  PoissonEdge1DSelfOrg::lambda0 = lambda0;
  PoissonEdge1DSelfOrg::lambdaMax= lambdaMax;
  // Initialise and run the sampler
  p_sampler = new smc::sampler<State>(numOfParticles, SMC_HISTORY_NONE);
  smc::moveset<State> Moveset(fInitialise, fMove, NULL);
  (*p_sampler).SetResampleParams(SMC_RESAMPLE_MULTINOMIAL, 1.01 * numOfParticles);
  (*p_sampler).SetMoveSet(Moveset);
  y[0] = initialObservation;
  (*p_sampler).Initialise();
}
