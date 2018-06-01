#include <R.h>
/* #include <apop.h> */
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include "functions.h"

/* LOOCV selects bandwidth for local linear kernel smoothing by leave-one-out cross validation. */

void LOOCV(double *Zin, int *nin, int *kin, double *LLKin, double *cv) {
  int const p=3;
  double const small_number=pow(10.0, -8);
  int i, j, i1, j1, n_local, ct;
  gsl_matrix *Z1, *wls_Xa, *wls_cova;
  gsl_vector *wls_w, *wls_ca, *wls_y;
  double leverage, r00, rssa, loocv, temp, x, y;

  /* Change the data type. */
  int n=nin[0], k=kin[0];
  gsl_matrix *Z, *LLK;

  Z = gsl_matrix_alloc(n, n);
  LLK = gsl_matrix_alloc(n, n);
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      gsl_matrix_set(Z, i, j, Zin[i+j*n]);
    }
  }
  

  /* Extend the observed image to avoid boundary problems. */
  Z1 = gsl_matrix_alloc(n+2*k, n+2*k);
  extend(Z, n, k, Z1);

  /* Compute common quantities. */
  i = k+3;
  j = k+3;
  n_local = 0;
  for (i1=(i-k); i1<=(i+k); i1++) {
    for (j1=(j-k); j1<=(j+k); j1++){
      if (((i1-i)*(i1-i) + (j1-j)*(j1-j)) <= k*k) {
  	n_local = n_local + 1;
      }
    }
  }

  wls_w = gsl_vector_alloc(n_local);
  wls_Xa = gsl_matrix_alloc(n_local, p);
  r00 = 0.0;
  ct = 0;
  for (i1=(i-k); i1<=(i+k); i1++){
    for (j1=(j-k); j1<=(j+k); j1++){
      if (((i1-i)*(i1-i) + (j1-j)*(j1-j)) <= k*k) {
  	x = (1.0*(i1-i))/n;
  	y = (1.0*(j1-j))/n;
  	gsl_matrix_set(wls_Xa, ct, 0, 1.0);
  	gsl_matrix_set(wls_Xa, ct, 1, x);
  	gsl_matrix_set(wls_Xa, ct, 2, y);
	temp = ker(x*n/k, y*n/k);
  	gsl_vector_set(wls_w, ct, temp);
	r00 = r00 + temp;
  	ct = ct + 1;
      }
    }
  }

  /* Compute Leverage */
  if (fabs(r00) <= small_number) {
    leverage = 1.0;
  } else {
    leverage = ker(0.0, 0.0)/r00;
  }

  if (fabs(leverage - 1.0) <= small_number) {
    error("leverage==1.0\n");
  }

  /* printf("n_local=%d, \t r00=%15.8f, \t leverage=%15.8f\n", n_local, r00, leverage); */

  /* Fit a local linear kernel regression. */
  loocv = 0.0;
  wls_y = gsl_vector_alloc(n_local);
  wls_ca = gsl_vector_alloc(p);
  wls_cova = gsl_matrix_alloc(p, p);
  /* LLK = gsl_matrix_alloc(n, n); */
  for (i=k; i<(n+k); i++) {
    for (j=k; j<(n+k); j++){
      ct = 0;
      for (i1=(i-k); i1<=(i+k); i1++) {
  	for (j1=(j-k); j1<=(j+k); j1++) {
  	  if (((i1-i)*(i1-i) + (j1-j)*(j1-j)) <= k*k) {
  	    gsl_vector_set(wls_y, ct, gsl_matrix_get(Z1, i1, j1));
  	    ct = ct + 1;
  	  }
  	}
      }
      if (ct != n_local) {
	error("ct != n_local\n");
      }
      gsl_multifit_linear_workspace *worka=gsl_multifit_linear_alloc(n_local,p);
      gsl_multifit_wlinear(wls_Xa, wls_w, wls_y, wls_ca, wls_cova, &rssa, worka);
      gsl_multifit_linear_free(worka);
      loocv = loocv + pow((gsl_vector_get(wls_ca, 0) - gsl_matrix_get(Z1, i, j))/(1.0 - leverage), 2);
      gsl_matrix_set(LLK, i-k, j-k, gsl_vector_get(wls_ca, 0));
    }
  }

  /* FILE *path_llk=fopen("LLK.out", "w"); */
  /* gsl_matrix_fprintf(path_llk, LLK, "%15.8f"); */
  /* fclose(path_llk); */
  
  loocv = loocv/n/n;

  /* Change the data type back. */
  cv[0] = loocv;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      LLKin[i+j*n] = gsl_matrix_get(LLK, i, j);
    }
  }
  gsl_matrix_free(Z);
  gsl_matrix_free(LLK);

  gsl_matrix_free(wls_cova);
  gsl_vector_free(wls_ca);
  gsl_vector_free(wls_y);
  gsl_vector_free(wls_w);
  gsl_matrix_free(wls_Xa);
  gsl_matrix_free(Z1);

}
