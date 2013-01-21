#ifndef ONLINEFITERING_H
#define ONLINEFITERING_H

#include "smctc.h"
#include "sampler.h"
#include "pfPoissonEdge1DSelfOrg.h"

using namespace std;
using namespace Rcpp;
using namespace PoissonEdge1DSelfOrg;

double onlineFiltering(smc::sampler<State>& Sampler, double observation, size_t i_obs);
void initialize(double initialObservation, int numOfParticles, double sigma_alpha, double alpha0, double beta0,
                            double lambda0, double lambdaMax);

#endif

