#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h> // included by R.h, so define USE_FC_LEN_T early

#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#ifndef FCONE
# define FCONE
#endif

void matsolve(double *a, double *b, int *nrowb, int *ncolb)
{
    int info;
    F77_CALL(dpotrf)("L", nrowb, a, nrowb, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed");
    F77_CALL(dpotrs)("L", nrowb, ncolb, a, nrowb, b, nrowb, &info FCONE);
    if (info != 0)
        error("solution failed");
}

void matdet(double *a, int *n, double *result)
{
    int info;
    F77_CALL(dpotrf)("L", n, a, n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed");
    int in = n[0];
    result[0] = 1.0;
    for (int i = 0; i < in; i++)
        result[0] *= a[i + in * i];
    result[0] *= result[0];
}
