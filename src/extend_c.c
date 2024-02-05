#include <R.h>

void extend_c(int *n_in, int *k_in, double *M, double *M1){
  int n = n_in[0];
  int k = k_in[0];
  int i, j;
  int n1 = n + 2*k;

  for (i = k; i < (k + n); i++){
    for (j = k; j < (k + n); j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, i-k, j-k)); */
      M1[i * n1 + j] = M[(i - k) * n + (j - k)];
    }
  }

  for (i = 0; i < k; i++){
    for (j = 0; j < k; j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, (k-1-j), (k-1-i))); */
      M1[i * n1 + j] = M[(k - 1 - j) * n + (k - 1 - i)];
    }
  }

  for (i = k; i < (k + n); i++){
    for (j = 0; j < k; j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, i-k, k-1-j)); */
      M1[i * n1 + j] = M[(i - k) * n + (k - 1 -j)];
    }
  }

  for (i = (k + n); i < (2 * k + n); i++){
    for (j = 0; j < k; j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, j+n-k, i-n-k)); */
      M1[i * n1 + j] = M[(j + n - k) * n + (i - n - k)];
    }
  }

  for (i = (k + n); i < (2 * k + n); i++){
    for (j = k; j < (k + n); j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, k+2*(n-1)+1-i, j-k)); */
      M1[i * n1 + j] = M[(k + 2 * (n - 1) + 1 - i) * n + (j - k)];
    }
  }

  for (i = (k + n); i < (2 * k + n); i++){
    for (j = (k + n); j < (2 * k + n); j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, k+2*(n-1)+1-j, k+2*(n-1)+1-i)); */
      M1[i * n1 + j] = M[(k + 2 * (n - 1) + 1 - j) * n + (k + 2 * (n - 1) + 1 - i)];
    }
  }

  for (i = k; i < (k + n); i++){
    for (j = (k + n); j < (2 * k + n); j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, i-k, k+2*(n-1)+1-j)); */
      M1[i * n1 + j] = M[(i - k) * n + (k + 2 * (n - 1) + 1 - j)];
    }
  }

  for (i = 0; i < k; i++){
    for (j = (k + n); j < (2 * k + n); j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, j-k-n, n-k+i)); */
      M1[i * n1 + j] = M[(j- k - n) * n + (n - k + i)];
    }
  }

  for (i = 0; i < k; i++){
    for (j = k; j < (k + n); j++){
      /* gsl_matrix_set(M1, i, j, gsl_matrix_get(M, k-1-i, j-k)); */
      M1[i * n1 + j] = M[(k - 1 - i) * n + (j - k)];
    }
  }
  
}
