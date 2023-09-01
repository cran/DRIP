#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Random.h>

/* call before generating any random variates */
void F77_SUB(rndstart)(void) { GetRNGstate(); }

/* call after done generating all random variates */
void F77_SUB(rndend)(void) { PutRNGstate(); }

/* call to generate one beta random variate */
double F77_SUB(myrunif)(double *alpha1, double *alpha2)
{
    return runif(alpha1[0], alpha2[0]);
}
