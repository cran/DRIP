#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "functions.h"

static R_NativePrimitiveArgType localdiff_types[4] = {INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType jp_llk_cv_types[5] = {INTSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType jp_llk_fit_types[6] = {INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType d_kq_types[4] = {INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType extend_types[4] = {INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType modify1_types[5] = {INTSXP, INTSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType modify2_types[4] = {INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType cluster_types[9] = {INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType cluster_denoise_bw_types[11] = {INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType cluster_deblur_bw_types[14] = {INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType roof_parsel_types[9] = {INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType roofdetect_types[6] = {INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType threestage_types[6] = {INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType threestage_denoise_bw_types[7] = {INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType threestage_deblur_bw_types[8] = {INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType jpex_loocv_types[5] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType JPEX0_types[7] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_FortranMethodDef fortranMethods[] = {
  {"lck_diff", (DL_FUNC) &F77_SUB(lck_diff), 4, localdiff_types},
  {"jp_llk_cv", (DL_FUNC) &F77_SUB(jp_llk_cv), 5, jp_llk_cv_types},
  {"jp_llk_fit", (DL_FUNC) &F77_SUB(jp_llk_fit), 6, jp_llk_fit_types},
  {"lc2k_diff", (DL_FUNC) &F77_SUB(lc2k_diff), 4, localdiff_types},
  {"d_kq", (DL_FUNC) &F77_SUB(d_kq), 4, d_kq_types},
  {"llk_diff", (DL_FUNC) &F77_SUB(llk_diff), 4, localdiff_types},
  {"ll2k_diff", (DL_FUNC) &F77_SUB(ll2k_diff), 4, localdiff_types},
  {"extend", (DL_FUNC) &F77_SUB(extend), 4, extend_types},
  {"modify1", (DL_FUNC) &F77_SUB(modify1), 5, modify1_types},
  {"modify2", (DL_FUNC) &F77_SUB(modify2), 4, modify2_types},
  {"cluster_cwm_denoise", (DL_FUNC) &F77_SUB(cluster_cwm_denoise), 9, cluster_types},
  {"cluster_cwm_deblur", (DL_FUNC) &F77_SUB(cluster_cwm_deblur), 9, cluster_types},
  {"cluster_cwm_denoise_bandwidth", (DL_FUNC) &F77_SUB(cluster_cwm_denoise_bandwidth), 11, cluster_denoise_bw_types},
  {"cluster_cwm_deblur_bandwidth", (DL_FUNC) &F77_SUB(cluster_cwm_deblur_bandwidth), 14, cluster_deblur_bw_types},
  {"roofdiff_denoise", (DL_FUNC) &F77_SUB(roofdiff_denoise), 4, localdiff_types},
  {"roofdiff_deblur", (DL_FUNC) &F77_SUB(roofdiff_deblur), 4, localdiff_types},
  {"roofedgeparsel_denoise", (DL_FUNC) &F77_SUB(roofedgeparsel_denoise), 9, roof_parsel_types},
  {"roofedgeparsel_deblur", (DL_FUNC) &F77_SUB(roofedgeparsel_deblur), 9, roof_parsel_types},
  {"roofdetect_denoise", (DL_FUNC) &F77_SUB(roofdetect_denoise), 6, roofdetect_types},
  {"roofdetect_deblur", (DL_FUNC) &F77_SUB(roofdetect_deblur), 6, roofdetect_types},
  {"deblur_3stage", (DL_FUNC) &F77_SUB(deblur_3stage), 6, threestage_types},
  {"denoise_3stage", (DL_FUNC) &F77_SUB(denoise_3stage), 6, threestage_types},
  {"denoise_3stage_bandwidth", (DL_FUNC) &F77_SUB(denoise_3stage_bandwidth), 7, threestage_denoise_bw_types},
  {"deblur_3stage_bandwidth", (DL_FUNC) &F77_SUB(deblur_3stage_bandwidth), 8, threestage_deblur_bw_types},
  {NULL, NULL, 0, NULL}
};

static R_CMethodDef cMethods[] = {
    {"LOOCV", (DL_FUNC) &LOOCV, 5, jpex_loocv_types},
    {"JPEX0", (DL_FUNC) &JPEX0, 7, JPEX0_types},
    {NULL, NULL, 0, NULL}
};

void attribute_visible R_init_DRIP(DllInfo* info) {
  R_registerRoutines(info, cMethods, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

