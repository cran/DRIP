#include <R.h>
#include <math.h>
#include <Rmath.h>
#include "functions.h"

/* LOOCV selects bandwidth for local linear kernel smoothing by leave-one-out cross validation. */

void LOOCV(double *Z, int *nin, int *kin, double *LLK, double *cv) {
  int p = 3, nc = 1;
  double const small_number = pow(10, -15);
  int i, j, i1, j1;
  double *Z1, *LLK_mat, *LLK_matdet, *LLK_coef;
  double A00, A10, A01, A11, A20, A02, dx, dy, temp, eta00, eta10, eta01, leverage, loocv;

  int n = nin[0], k = kin[0];
  double h = (1.0 * k)/n;
 

  /* Extend the observed image to avoid boundary problems. */

  Z1 = (double *)malloc((n + 2 * k) * (n + 2 * k) * sizeof(double));
  extend_c(nin, kin, Z, Z1);

  /* Compute common quantities. */

  i = k+3;
  j = k+3;
  A00 = 0.0;
  A10 = 0.0;
  A01 = 0.0;
  A11 = 0.0;
  A20 = 0.0;
  A02 = 0.0;

  for (i1 = (i - k); i1 <= (i + k); i1++) {
    for (j1 = (j - k); j1 <= (j + k); j1++){

      if (((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) < (k * k)) {
	dy = -(i1 - i)/(1.0 * n);/* recall the way the matrix is entered and the way the image is expressed in terms of x and y */
	dx = (j1 - j)/(1.0 * n);
	temp = ker(dx/h, dy/h);
	A00 = A00 + temp;
	A10 = A10 + dx * temp;
	A01 = A01 + dy * temp;
	A11 = A11 + dx * dy * temp;
	A20 = A20 + dx * dx * temp;
	A02 = A02 + dy * dy * temp;
      }

      
    }
  }

  LLK_mat = (double *)malloc(p * p * sizeof(double));
  LLK_mat[0] = A00;
  LLK_mat[1] = A10;
  LLK_mat[2] = A01;
  LLK_mat[3] = A10;
  LLK_mat[4] = A20;
  LLK_mat[5] = A11;
  LLK_mat[6] = A01;
  LLK_mat[7] = A11;
  LLK_mat[8] = A02;

  if (fabs(A00) < small_number) {
    error("The bandwidth is too small in LOOCV. \n");
  }

  LLK_matdet = (double *)malloc(1 * sizeof(double));
  matdet(LLK_mat, &p, LLK_matdet);

  if(LLK_matdet[0] < small_number) {
    error("The bandwidth is too small for LLK smoothing in LOOCV. \n");
  }

  /* Compute Leverage */
  
  leverage = ker(0.0, 0.0)/A00;


  /* Fit a local linear kernel regression. */

  loocv = 0.0; 
  LLK_coef = (double *)malloc(p * sizeof(double));

  for (i = k; i < (n + k); i++) {
    for (j = k; j < (n + k); j++){

      eta00 = 0.0;
      eta10 = 0.0;
      eta01 = 0.0;

      for (i1 = (i - k); i1 <= (i + k); i1++) {
	for (j1 = (j - k); j1 <= (j + k); j1++) {
	  if (((i1 - i) * (i1 - i) + (j1 - j) *(j1 - j)) < (k * k)) {
	    dy = -(i1 - i)/(1.0 * n); /* recall the way the matrix is entered and the way the image is expressed in terms of x and y */
	    dx = (j1 - j)/(1.0 * n);
	    temp = Z1[i1 * (n + 2 * k) + j1] * ker(dx/h, dy/h);
	    eta00 = eta00 + temp;
	    eta10 = eta10 + dx * temp;
	    eta01 = eta01 + dy * temp;

	  }
	}
      }
      /* Need to reset the values in LLK_mat because it gets changed everytime we run matsolve */
      LLK_mat[0] = A00;
      LLK_mat[1] = A10;
      LLK_mat[2] = A01;
      LLK_mat[3] = A10;
      LLK_mat[4] = A20;
      LLK_mat[5] = A11;
      LLK_mat[6] = A01;
      LLK_mat[7] = A11;
      LLK_mat[8] = A02;

      

      LLK_coef[0] = eta00;
      LLK_coef[1] = eta10;
      LLK_coef[2] = eta01;

      matsolve(LLK_mat, LLK_coef, &p, &nc);

      LLK[(i - k) * n + (j - k)] = LLK_coef[0];

      loocv = loocv + pow((LLK_coef[0] - Z1[i * (n + 2 * k) + j])/(1.0 - leverage), 2);
      
    }
  }

  loocv = loocv/n/n;
  cv[0] = loocv;

  free(Z1);
  free(LLK_mat);
  free(LLK_matdet);
  free(LLK_coef);

}
