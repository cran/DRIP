#ifndef DRIP_FUNCTIONS_H
#define DRIP_FUNCTIONS_H

#include <R.h>
#include <Rinternals.h>

void F77_SUB(lck_diff)(int *n_in, double *obsImg_in, int *bandwidth_in, double *diff_in);
void F77_SUB(jp_llk_cv)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidth_in, double *cv_in);
void F77_SUB(jp_llk_fit)(int *n_in, double *obsImg_in, int *bandwidth_in, double *fitted_in,
			 double *resid_in, double *sigma_in);
void F77_SUB(lc2k_diff)(int *n_in, double *obsImg_in, int *bandwidth_in, double *diff_in);
void F77_SUB(d_kq)(int *n_in, int *edge1_in, int *edge2_in, double *dkq_in);
void F77_SUB(llk_diff)(int *n_in, double *obsImg_in, int *bandwidth_in, double *diff_in);
void F77_SUB(ll2k_diff)(int *n_in, double *obsImg_in, int *bandwidth_in, double *diff_in);
void F77_SUB(extend)(int *n_in, int *k_in, double *z_in, double *z1_in);
void F77_SUB(modify1)(int *n_in, int *k_in, int *bound_in, double *z_in, int *edge_in);
void F77_SUB(modify2)(int *n_in, int *k_in, int *bound_in, int *edge_in);
void F77_SUB(cluster_cwm_denoise)(int *n_in, double *obsImg_in, int *k_in, double *zq_in, double *sigma_in,
				  double *phi0_in, double *mean_std_abs_in, int *cw_in, double *estImg_in);
void F77_SUB(cluster_cwm_deblur)(int *n_in, double *obsImg_in, int *k_in, double *zq_in, double *sigma_in,
				 double *phi0_in, double *mean_std_abs_in, int *cw_in, double *estImg_in);
void F77_SUB(cluster_cwm_denoise_bandwidth)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidths_in,
					    double *zq_in, double *sigma_in, double *phi0_in,
					    double *mean_std_abs_in, int *cw_in, int *bandwidth_hat_in,
					    double *cv_in);
void F77_SUB(cluster_cwm_deblur_bandwidth)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidths_in,
					   double *zq_in, double *sigma_in, double *phi0_in,
					   double *mean_std_abs_in, int *cw_in, double *relwt_in,
					   int *bandwidth_hat_in, double *cv_in, double *cv_jump_in,
					   double *cv_cty_in);
void F77_SUB(roofdiff_denoise)(int *n_in, double *obsImg_in, int *bandwidth_in, double *diff_in);
void F77_SUB(roofdiff_deblur)(int *n_in, double *obsImg_in, int *bandwidth_in, double *diff_in);
void F77_SUB(roofedgeparsel_denoise)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidth_in,
				     int *nthresh_in, double *thresh_in, int *nboot_in, int *edge1_in,
				     double *dKQ_in);
void F77_SUB(roofedgeparsel_deblur)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidth_in,
				     int *nthresh_in, double *thresh_in, int *nboot_in, int *edge1_in,
				     double *dKQ_in);
void F77_SUB(roofdetect_denoise)(int *n_in, double *obsImg_in, int *bandwidth_in, double *thresh_in,
				 int *edge1_in, int *edge2_in);
void F77_SUB(roofdetect_deblur)(int *n_in, double *obsImg_in, int *bandwidth_in, double *thresh_in,
				 int *edge1_in, int *edge2_in);
void F77_SUB(deblur_3stage)(int *n_in, double *obsImg_in, int *bandwidth_in, int *edge1_in, int *edge2_in,
			    double *estImg_in);
void F77_SUB(denoise_3stage)(int *n_in, double *obsImg_in, int *bandwidth_in, int *edge1_in, int *edge2_in,
			    double *estImg_in);
void F77_SUB(denoise_3stage_bandwidth)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidth_in,
				       int *edge1_in, int *edge2_in, double *cv_in);
void F77_SUB(deblur_3stage_bandwidth)(int *n_in, double *obsImg_in, int *nband_in, int *bandwidth_in,
				      int *edge1_in, int *edge2_in, int *nboot_in, double *msecv_in);


void LOOCV(double *Z, int *nin, int *kin, double *LLK, double *cv);
void JPEX0(double *Z, int *nin, int *kin, double *alphain, double *sigmain, double *EDGE, double *fhat);  

extern void extend_c(int *n_in, int *k_in, double *M, double *M1);
extern double ker(double x, double y);
extern void matsolve(double *a, double *b, int *nrowb, int *ncolb);
extern void matdet(double *a, int *n, double *result);

#endif /* DRIP_FUNCTIONS_H */
