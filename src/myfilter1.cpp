#include <Rcpp.h>
#include "smctc.h"
#include <iostream>
#include <cmath>
//#include <gsl/gsl_randist.h>
using namespace std;
using namespace Rcpp;

#include "myfilter1.h"

// global variables
std::vector<cv_obs> yy;

double integrand_mean_x(const cv_state&, void*);
double integrand_mean_sigma_x(const cv_state&, void*);
double integrand_mean_sigma_y(const cv_state&, void*);
double integrand_mean_y(const cv_state&, void*);
double integrand_var_x(const cv_state&, void*);
double integrand_var_y(const cv_state&, void*);

// [[Rcpp::export]]
SEXP myfilter1(unsigned long particleNum, NumericVector observation) { 	

    long lIterates;

    try {
        lIterates = observation.size();
        yy.reserve(lIterates);
        for (unsigned long i = 0; i < lIterates; ++i) {
          (yy[i]).val = (double)(observation[i]);
        }

        //Initialise and run the sampler
        smc::sampler<cv_state> Sampler(particleNum, SMC_HISTORY_NONE);  
        smc::moveset<cv_state> Moveset(fInitialise, fMove, NULL);

        Sampler.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();

        Rcpp::NumericVector Xm(lIterates), Xv(lIterates), sigma_x_m(lIterates), sigma_y_m(lIterates);

        for(int n=0; n < lIterates; ++n) {
            Sampler.Iterate();
            Xm(n) = Sampler.Integrate(integrand_mean_x, NULL);
            Xv(n) = Sampler.Integrate(integrand_var_x, (void*)&Xm(n));
            sigma_x_m(n) = Sampler.Integrate(integrand_mean_sigma_x, NULL);
            sigma_y_m(n) = Sampler.Integrate(integrand_mean_sigma_y, NULL);
        }

        return Rcpp::DataFrame::create(Rcpp::Named("Xm") = Xm,
                                       Rcpp::Named("Xv") = Xv,
                                       Rcpp::Named("sxm") = sigma_x_m,
                                       Rcpp::Named("sym") = sigma_y_m);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o 
        //exit(e.lCode);		// we're just called from R so we should not exit
    }
    return R_NilValue;          	// to provide a return 
}

double integrand_mean_x(const cv_state& s, void *)
{
  return s.val;
}
double integrand_mean_sigma_x(const cv_state& s, void *)
{
  return s.x_sigma;
}
double integrand_mean_sigma_y(const cv_state& s, void *)
{
  return s.y_sigma;
}
double integrand_var_x(const cv_state& s, void* vmx)
{
  double* dmx = (double*)vmx;
  double d = (s.val - (*dmx));
  return d*d;
}

const double var_s0 = 4;
double var_s  = 0.1;
const double var_state_x_sigma0 = 4;
const double var_state_x_sigma = 0.5;
double var_state_y_sigma0 = 4;
double var_state_y_sigma = 0.5;
///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
double logLikelihood(long lTime, const cv_state & X)
{
    return -0.5 * pow((yy[lTime].val - X.val),2) / X.y_sigma;
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<cv_state> fInitialise(smc::rng *pRng)
{
  cv_state value;
  value.val = pRng->Normal(0,sqrt(var_s0));
  value.x_sigma = sqrt(var_state_x_sigma0);
  value.y_sigma = sqrt(var_state_y_sigma0);
  return smc::particle<cv_state>(value,logLikelihood(0,value));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng)
{
  try{
  cv_state * cv_to = pFrom.GetValuePointer();
  cv_to->val += pRng->Normal(0,cv_to->x_sigma);
  cv_to->x_sigma += pRng->Normal(0,sqrt(var_state_x_sigma));
  cv_to->y_sigma += pRng->Normal(0,sqrt(var_state_y_sigma));
  pFrom.AddToLogWeight(logLikelihood(lTime, *cv_to));
  } catch (smc::exception e){
    
  }
}
