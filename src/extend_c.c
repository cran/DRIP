#include <gsl/gsl_matrix.h>

void extend(gsl_matrix *M, int const n, int const k, gsl_matrix *M1){
  int i, j;

  for (i=k; i<(k+n); i++){
    for (j=k; j<(k+n); j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, i-k, j-k));
    }
  }

  for (i=0; i<k; i++){
    for (j=0; j<k; j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, (k-1-j), (k-1-i)));
    }
  }

  for (i=k; i<(k+n); i++){
    for (j=0; j<k; j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, i-k, k-1-j));
    }
  }

  for (i=(k+n); i<(2*k+n); i++){
    for (j=0; j<k; j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, j+n-k, i-n-k));
    }
  }

  for (i=(k+n); i<(2*k+n); i++){
    for (j=k; j<(k+n); j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, k+2*(n-1)+1-i, j-k));
    }
  }

  for (i=(k+n); i<(2*k+n); i++){
    for (j=(k+n); j<(2*k+n); j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, k+2*(n-1)+1-j, k+2*(n-1)+1-i));
    }
  }

  for (i=k; i<(k+n); i++){
    for (j=(k+n); j<(2*k+n); j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, i-k, k+2*(n-1)+1-j));
    }
  }

  for (i=0; i<k; i++){
    for (j=(k+n); j<(2*k+n); j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, j-k-n, n-k+i));
    }
  }

  for (i=0; i<k; i++){
    for (j=k; j<(k+n); j++){
      gsl_matrix_set(M1, i, j, gsl_matrix_get(M, k-1-i, j-k));
    }
  }
  
}
