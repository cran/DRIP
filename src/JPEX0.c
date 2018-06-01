#include <R.h>
/* #include <apop.h> */
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include "functions.h"

/* JPEX0 deblurs the image using jump-preserving constant extrapolation. */

void JPEX0(double *Zin, int *nin, int *kin, double *alphain, double *sigmain, double *EDGEin, double *fhatin){
  int const p=3, p0=1;
  int const rex=25; 		/* radius for searching the nearest sharp pixel */
  int df, itemp, jtemp;
  int i, j, i1, j1, n_local, ct;
  size_t istar, jstar;
  double x, y, rss0, rssa, thresh;
  gsl_matrix *Z1, *Ttilde, *wls_cov0, *wls_cova, *wls_X0, *wls_Xa, *dist;
  gsl_vector *wls_w, *wls_c0, *wls_ca, *wls_y;

  /* Change the data types */
  int n=nin[0], k=kin[0];
  double alpha=alphain[0], sigma=sigmain[0];
  double sigma2=sigma*sigma;
  gsl_matrix *Z, *EDGE, *fhat;

  Z = gsl_matrix_alloc(n, n);
  EDGE = gsl_matrix_alloc(n, n);
  fhat = gsl_matrix_alloc(n, n);
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      gsl_matrix_set(Z, i, j, Zin[i+j*n]);
      gsl_matrix_set(EDGE, i, j, EDGEin[i+j*n]);
      gsl_matrix_set(fhat, i, j, fhatin[i+j*n]);
    }
  }
  
  
  /* Extend the observed image to avoid boundary problems. */
  Z1 = gsl_matrix_alloc(n+2*k, n+2*k);
  extend(Z, n, k, Z1);

  /* Flag blurry pixels. */

  /* Compute common quantities. */
  i = k+3;
  j = k+3;
  n_local = 0;
  for (i1=(i-k); i1<=(i+k); i1++) {
    for (j1=(j-k); j1<=(j+k); j1++){
      if (((1.0*i1-1.0*i)*(1.0*i1-1.0*i) + (1.0*j1-1.0*j)*(1.0*j1-1.0*j)) <= k*k) {
  	n_local = n_local + 1;
      }
    }
  }

  wls_w = gsl_vector_alloc(n_local);
  wls_X0 = gsl_matrix_alloc(n_local, p0);
  wls_Xa = gsl_matrix_alloc(n_local, p);
  ct = 0;
  for (i1=(i-k); i1<=(i+k); i1++){
    for (j1=(j-k); j1<=(j+k); j1++){
      if (((1.0*i1-1.0*i)*(1.0*i1-1.0*i) + (1.0*j1-1.0*j)*(1.0*j1-1.0*j)) <= k*k) {
  	x = (1.0*i1-1.0*i)/n;
  	y = (1.0*j1-1.0*j)/n;
	gsl_matrix_set(wls_X0, ct, 0, 1.0);
  	gsl_matrix_set(wls_Xa, ct, 0, 1.0);
  	gsl_matrix_set(wls_Xa, ct, 1, x);
  	gsl_matrix_set(wls_Xa, ct, 2, y);
  	gsl_vector_set(wls_w, ct, ker(x*n/k, y*n/k));
  	ct = ct + 1;
      }
    }
  }
  /* printf("Xa[20,0]=%15.8f\n", gsl_matrix_get(wls_Xa, 20, 0)); */
  /* printf("Xa[20,1]=%15.8f\n", gsl_matrix_get(wls_Xa, 20, 1)); */
  /* printf("Xa[20,2]=%15.8f\n", gsl_matrix_get(wls_Xa, 20, 2)); */
  /* printf("w[20]=%15.8f\n", gsl_vector_get(wls_w, 20)); */
  wls_y = gsl_vector_alloc(n_local);
  wls_c0 = gsl_vector_alloc(p0);
  wls_ca = gsl_vector_alloc(p);
  wls_cov0 = gsl_matrix_alloc(p0,p0);
  wls_cova = gsl_matrix_alloc(p, p);

  /* Perform Chi square test for each neighborhood. */
  df = p - p0;
  thresh = gsl_cdf_chisq_Pinv(1.0-alpha, df);
  Ttilde = gsl_matrix_alloc(n, n);
  for (i=k; i<(n+k); i++) {
    for (j=k; j<(n+k); j++){
      ct = 0;
      for (i1=0; i1<=(2*k); i1++) {
  	for (j1=0; j1<=(2*k); j1++) {
	  itemp = i + k - (2*k - i1);
	  jtemp = j + k - (2*k - j1);
  	  if (((itemp-1.0*i)*(itemp-1.0*i) + (jtemp-1.0*j)*(jtemp-1.0*j)) <= k*k) {
  	    gsl_vector_set(wls_y, ct, gsl_matrix_get(Z1, itemp, jtemp));
  	    ct = ct + 1;
  	  }
  	}
      }
      gsl_multifit_linear_workspace *work0=gsl_multifit_linear_alloc(n_local,p0);
      gsl_multifit_wlinear(wls_X0, wls_w, wls_y, wls_c0, wls_cov0, &rss0, work0);
      itemp = i - k;
      jtemp = j - k;
      gsl_matrix_set(Ttilde, itemp, jtemp, gsl_vector_get(wls_c0, 0));
      gsl_multifit_linear_free(work0);
      gsl_multifit_linear_workspace *worka=gsl_multifit_linear_alloc(n_local,p);
      gsl_multifit_wlinear(wls_Xa, wls_w, wls_y, wls_ca, wls_cova, &rssa, worka);
      gsl_multifit_linear_free(worka);
      gsl_matrix_set(EDGE, itemp, jtemp, (rss0 - rssa)/(sigma2));
      /* if (((rss0 - rssa)/(sigma*sigma)) >= thresh) { */
      /* 	gsl_matrix_set(EDGE, i-k, j-k, 1); */
      /* } else { */
      /* 	gsl_matrix_set(EDGE, i-k, j-k, 0); */
      /* } */
    }
  }

  /* printf("Tilde[0,0]=%15.8f\n", gsl_matrix_get(Ttilde, 0, 0)); */
  /* printf("Tilde[20,20]=%15.8f\n", gsl_matrix_get(Ttilde, 20, 0)); */
  /* printf("Tilde[254,254]=%15.8f\n", gsl_matrix_get(Ttilde, 254, 254)); */
  /* printf("wls_y[40]=%15.8f\n", gsl_vector_get(wls_y, 40)); */


  /* Extrapolate from sharp pixels. */
  /* gsl_matrix *istarmap1, *jstarmap1, *istarmap2, *jstarmap2; */
  /* istarmap1 = gsl_matrix_alloc(n, n); */
  /* jstarmap1 = gsl_matrix_alloc(n, n); */
  /* istarmap2 = gsl_matrix_alloc(n, n); */
  /* jstarmap2 = gsl_matrix_alloc(n, n); */

  /* dist = gsl_matrix_alloc(n, n); */
  /* for (i=0; i<n; i++) { */
  /*   for (j=0; j<n; j++) { */
  /*     gsl_matrix_set_all(dist, INFINITY); */
  /*     if (gsl_matrix_get(EDGE, i, j) >= thresh) { */
  /*       for (i1=0; i1<n; i1++) { */
  /*         for (j1=0; j1<n; j1++) { */
  /*           if (gsl_matrix_get(EDGE, i1, j1) < thresh) { */
  /*             gsl_matrix_set(dist, i1, j1, 1.0*(i1-i)*(i1-i) + (j1-j)*(j1-j)); */
  /*           } */
  /*         } */
  /*       } */
  /* 	gsl_matrix_min_index(dist, &istar, &jstar); */
  /* 	if (i==0 && j==0) { */
  /* 	  printf("istar=%zu, \t jstar=%zu\n", istar, jstar); */
  /* 	} */
  /*       gsl_matrix_set(istarmap1, i, j, istar); */
  /* 	gsl_matrix_set(jstarmap1, i, j, jstar); */
  /*       gsl_matrix_set(fhat, i, j, gsl_matrix_get(Ttilde, istar, jstar)); */
  /*     } else { */
  /*       gsl_matrix_set(fhat, i, j, gsl_matrix_get(Ttilde, i, j)); */
  /* 	gsl_matrix_set(istarmap1, i, j, i); */
  /* 	gsl_matrix_set(jstarmap1, i, j, j); */
  /*     } */
  /*   } */
  /* } */

  /* gsl_matrix_free(dist); */
  dist = gsl_matrix_alloc(2*rex+1, 2*rex+1);
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      gsl_matrix_set_all(dist, INFINITY);
      if (gsl_matrix_get(EDGE, i, j) >= thresh) {
  	for (i1=0; i1<(2*rex+1); i1++) {
  	  for (j1=0; j1<(2*rex+1); j1++) {
  	    if (i+rex-i1<n && i+rex-i1>=0 && j+rex-j1>=0 && j+rex-j1<n) {
  	      if (gsl_matrix_get(EDGE, i+rex-i1, j+rex-j1) < thresh) {
  		gsl_matrix_set(dist, i1, j1, 1.0*(rex-i1)*(rex-i1) + (rex-j1)*(rex-j1));
  	      }
  	    }
  	  }
  	}
  	gsl_matrix_min_index(dist, &istar, &jstar);
	/* if (i==0 && j==0) { */
	/*   printf("istar=%d, \t jstar=%d\n", i+rex-((int)istar), j+rex-((int)jstar)); */
	/* } */
	/* gsl_matrix_set(istarmap2, i, j, i+rex-((int)istar)); */
	/* gsl_matrix_set(jstarmap2, i, j, j+rex-((int)jstar)); */
  	gsl_matrix_set(fhat, i, j, gsl_matrix_get(Ttilde, i+rex-((int)istar), j+rex-((int)jstar)));
      } else {
  	gsl_matrix_set(fhat, i, j, gsl_matrix_get(Ttilde, i, j));
	/* gsl_matrix_set(istarmap2, i, j, i); */
	/* gsl_matrix_set(jstarmap2, i, j, j); */
      }
    }
  }

  /* FILE *istar1=fopen("istar1.out", "w"); */
  /* gsl_matrix_fprintf(istar1, istarmap1, "%15.8f"); */
  /* fclose(istar1); */
  /* FILE *jstar1=fopen("jstar1.out", "w"); */
  /* gsl_matrix_fprintf(jstar1, jstarmap1, "%15.8f"); */
  /* fclose(jstar1); */
  /* FILE *istar2=fopen("istar2.out", "w"); */
  /* gsl_matrix_fprintf(istar2, istarmap2, "%15.8f"); */
  /* fclose(istar2); */
  /* FILE *jstar2=fopen("jstar2.out", "w"); */
  /* gsl_matrix_fprintf(jstar2, jstarmap2, "%15.8f"); */
  /* fclose(jstar2); */

  /* gsl_matrix_free(istarmap1); */
  /* gsl_matrix_free(jstarmap1); */
  /* gsl_matrix_free(istarmap2); */
  /* gsl_matrix_free(jstarmap2); */

  /* Change the data type back. */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      fhatin[i+j*n] = gsl_matrix_get(fhat, i, j);
      EDGEin[i+j*n] = gsl_matrix_get(EDGE, i, j);
    }
  }
  gsl_matrix_free(Z);
  gsl_matrix_free(EDGE);
  gsl_matrix_free(fhat);

  gsl_matrix_free(wls_cova);
  gsl_matrix_free(wls_Xa);
  gsl_vector_free(wls_ca);
  gsl_matrix_free(wls_cov0);
  gsl_matrix_free(wls_X0);
  gsl_vector_free(wls_c0);
  gsl_vector_free(wls_y);
  gsl_vector_free(wls_w);
  gsl_matrix_free(Ttilde);
  gsl_matrix_free(dist);

}
