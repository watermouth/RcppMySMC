#include "smctc.h"

class cv_state 
{
public:
    double val;
};

class cv_obs
{
public:
    double val;
};

double logLikelihood(long lTime, const cv_state & X);

smc::particle<cv_state> fInitialise(smc::rng *pRng);
void fMove(long lTime, smc::particle<cv_state> & pFrom, smc::rng *pRng);

extern std::vector<cv_obs> yy;
