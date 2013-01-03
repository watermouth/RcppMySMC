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
double integrand_mean_y(const cv_state&, void*);
double integrand_var_x(const cv_state&, void*);
double integrand_var_y(const cv_state&, void*);

// [[Rcpp::export]]
SEXP myfilter1(unsigned long particleNum, NumericVector observation, bool useF, Function f) { 	

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

        Rcpp::NumericVector Xm(lIterates), Xv(lIterates), Ym(lIterates), Yv(lIterates);

        for(int n=0; n < lIterates; ++n) {
            Sampler.Iterate();
      
            Xm(n) = Sampler.Integrate(integrand_mean_x, NULL);
            Xv(n) = Sampler.Integrate(integrand_var_x, (void*)&Xm(n));
            Ym(n) = Sampler.Integrate(integrand_mean_y, NULL);
            Yv(n) = Sampler.Integrate(integrand_var_y, (void*)&Ym(n));

            if (useF) f(Xm, Ym);
        }

        return Rcpp::DataFrame::create(Rcpp::Named("Xm") = Xm,
                                       Rcpp::Named("Xv") = Xv,
                                       Rcpp::Named("Ym") = Ym,
                                       Rcpp::Named("Yv") = Yv);
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

double integrand_var_x(const cv_state& s, void* vmx)
{
  double* dmx = (double*)vmx;
  double d = (s.val - (*dmx));
  return d*d;
}

double integrand_mean_y(const cv_state& s, void *)
{
  return s.val;
}

double integrand_var_y(const cv_state& s, void* vmy)
{
  double* dmy = (double*)vmy;
  double d = (s.val - (*dmy));
  return d*d;
}

double var_s0 = 4;
double var_u0 = 1;
double var_s  = 0.1;
double var_u  = 0.1;

///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
double logLikelihood(long lTime, const cv_state & X)
{
    return -0.5 * pow((yy[lTime].val - X.val),2) / var_u;
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<cv_state> fInitialise(smc::rng *pRng)
{
  cv_state value;
  value.val = pRng->Normal(0,sqrt(var_s0));
  return smc::particle<cv_state>(value,logLikelihood(0,value));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng)
{
  cv_state * cv_to = pFrom.GetValuePointer();
  cv_to->val += pRng->Normal(0,sqrt(var_s));
  pFrom.AddToLogWeight(logLikelihood(lTime, *cv_to));
}
