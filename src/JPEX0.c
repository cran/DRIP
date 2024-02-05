#include <R.h>
#include <math.h>
#include <Rmath.h>
#include "functions.h"

int minloc(int n, double *x) {
  int i, min_loc;
  double xmin;

  xmin = x[0];
  min_loc = 0;
  for (i = 0; i < n; i++) {
    if (xmin > x[i]) {
      xmin = x[i];
      min_loc = i;
    }
  }
  return(min_loc);
}

/* JPEX0 deblurs the image using jump-preserving constant extrapolation. */

void JPEX0(double *Z, int *nin, int *kin, double *alphain, double *sigmain, double *EDGE, double *fhat){
  int p = 3, nc = 1, p0 = 1;
  /* int rex = 25; */ 		/* radius for searching the nearest sharp pixel */
  int df, i, j, i1, j1;
  int starLoc, istar, jstar;
  double rss0, rssa, thresh;
  double *Z1, *Ttilde, *dist, *LCK_coef, *LLK_coef, *LLK_mat, *LLK_matdet;
  double A00, A10, A01, A11, A20, A02, dx, dy, temp, eta00, eta10, eta01, S00;

  double const small_number = pow(10, -15);

  /* Change the data types */
  int n = nin[0], k = kin[0];
  double alpha = alphain[0], sigma = sigmain[0];
  double sigma2 = sigma * sigma;
  double h = (1.0 * k)/n;

  /* Extend the observed image to avoid boundary problems. */

  Z1 = (double *)malloc((n + 2 * k) * (n + 2 * k) * sizeof(double));
  extend_c(nin, kin, Z, Z1);

  /* Flag blurry pixels. */

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

      if (((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= (k * k)) {
	dy = (i1 - i)/(1.0 * n); /* Recall the way a matrix from R is stacked in C and the way an image is expressed in terms of x and y. */
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

  /* Rprintf("%f, %f, %f \n", A00, A10, A01); */
  /* Rprintf("%f, %f, %f \n", A10, A20, A11); */
  /* Rprintf("%f, %f, %f \n", A01, A11, A02); */

  if (fabs(A00) < small_number) {
    error("The bandwidth is too small for LCK smoothing. \n");
  }

  LLK_matdet = (double *)malloc(1 * sizeof(double));
  matdet(LLK_mat, &p, LLK_matdet);

  if(LLK_matdet[0] < small_number) {
    error("The bandwidth is too small for LLK smoothing. \n");
  }
  

  LCK_coef = (double *)malloc(p0 * sizeof(double));
  LLK_coef = (double *)malloc(p * sizeof(double));

  /* Perform Chi square test for each neighborhood. */
  df = p - p0;
  thresh = qchisq(1.0 - alpha, 1.0 * df, 1, 0);
  Ttilde = (double *)malloc(n * n * sizeof(double));

  for (i = k; i < (n + k); i++) {
    for (j = k; j < (n + k); j++){

      eta00 = 0.0;
      eta10 = 0.0;
      eta01 = 0.0;
      S00 = 0.0;
      for (i1 = (i - k); i1 <= (i + k); i1++) {
	for (j1 = (j - k); j1 <= (j + k); j1++) {
	  if (((i1 - i) * (i1 - i) + (j1 - j) *(j1 - j)) <= (k * k)) {
	    dy = (i1 - i)/(1.0 * n); /* Recall the way a matrix from R is stacked in C and the way an image is expressed in terms of x and y. */
	    dx = (j1 - j)/(1.0 * n);
	    temp = Z1[i1 * (n + 2 * k) + j1] * ker(dx/h, dy/h);
	    S00 = S00 + Z1[i1 * (n + 2 * k) + j1] * temp;
	    eta00 = eta00 + temp;
	    eta10 = eta10 + dx * temp;
	    eta01 = eta01 + dy * temp;
	  }
	}
      }
      LCK_coef[0] = eta00/A00;
      rss0 = S00 - 2 * LCK_coef[0] * eta00 + pow(LCK_coef[0], 2) * A00;

      LLK_coef[0] = eta00;
      LLK_coef[1] = eta10;
      LLK_coef[2] = eta01;
      /* Need to reset values in LLK_mat because of the way it works. */
      LLK_mat[0] = A00;
      LLK_mat[1] = A10;
      LLK_mat[2] = A01;
      LLK_mat[3] = A10;
      LLK_mat[4] = A20;
      LLK_mat[5] = A11;
      LLK_mat[6] = A01;
      LLK_mat[7] = A11;
      LLK_mat[8] = A02;
      matsolve(LLK_mat, LLK_coef, &p, &nc);
      rssa = S00 + LLK_coef[0] * LLK_coef[0] * A00 + LLK_coef[1] * LLK_coef[1] * A20 + LLK_coef[2] * LLK_coef[2] * A02 \
	- 2 * (LLK_coef[0] * eta00 + LLK_coef[1] * eta10 + LLK_coef[2] * eta01) \
	+ 2 * (LLK_coef[0] * LLK_coef[1] * A10 + LLK_coef[0] * LLK_coef[2] * A01 + LLK_coef[1] * LLK_coef[2] * A11);

      EDGE[(i - k) * n + (j - k)] = (rss0 - rssa)/sigma2;
      Ttilde[(i - k) * n + (j - k)] = LCK_coef[0];

    }
  }


  /* Extrapolate from sharp pixels. */

  /* dist = (double *)malloc((2 * rex + 1) * (2 * rex + 1) * sizeof(double)); */
  
  dist = (double *)malloc(n * n * sizeof(double));
  
  
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      
      /* for (i1 = 0; i1 < ((2 * rex + 1) * (2 * rex + 1)); i1++) { */
      /* 	dist[i1] = INFINITY; */
      /* } */
      /* edge_stat = EDGE[i * n + j]; */
      
      if (EDGE[i * n + j] >= thresh) {

	for (i1 = 0; i1 < (n * n); i1++) {
	  dist[i1] = INFINITY;
	}

  	for (i1 = 0; i1 < n; i1++) {
  	  for (j1 = 0; j1 < n; j1++) {

	    if (EDGE[i1 * n + j1] < thresh) {

	      dist[i1 * n + j1] = 1.0 * (i - i1) * (i - i1) + (j - j1) * (j - j1);

	    }

	    /* if (i + rex - i1 < n && i + rex - i1 >= 0 && j + rex - j1>= 0 && j + rex -j1 < n) { */

	    /*   itemp = i + rex - i1; */
	    /*   jtemp = j + rex - j1; */
  	    /*   if (EDGE[itemp * n + jtemp] < thresh) { */

	    /* 	/\* gsl_matrix_set(dist, i1, j1, 1.0*(rex-i1)*(rex-i1) + (rex-j1)*(rex-j1)); *\/ */
	    /* 	dist[i1 * (2 * rex + 1) + j1] = 1.0 * (rex - i1) * (rex - i1) + (rex - j1) * (rex - j1); */

	    /*   } */
  	    /* } */
  	  }
  	}

	/* starLoc = minloc((2 * rex + 1) * (2 * rex + 1), dist); */
	/* istar = starLoc / (2 * rex + 1); */
	/* jstar = starLoc % (2 * rex + 1); */
	/* fhat[i * n + j] = Ttilde[(i + rex - istar) * n + (j + rex - jstar)]; */

	starLoc = minloc(n * n, dist);
	istar = starLoc / n;
	jstar = starLoc % n;
	fhat[i * n + j] = Ttilde[istar * n + jstar];
	
      } else {

	fhat[i * n + j] = Ttilde[i * n + j];

      }
      
    }
  }

 
  free(LLK_matdet);
  free(LLK_mat);
  free(LLK_coef);
  free(LCK_coef);
  free(Z1);
  free(Ttilde);
  free(dist);

}
