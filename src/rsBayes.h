#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP phenoMod0(SEXP y_r, SEXP t_r, SEXP n_r, SEXP family_r, SEXP tNormBounds_r,
		 SEXP gamma_r,
                 SEXP phiIG_r, SEXP alphaUa_r, SEXP alphaUb_r,
                 SEXP phiStarting_r, SEXP alphaStarting_r,
                 SEXP phiTuning_r, SEXP alphaTuning_r,
                 SEXP nSamples_r, SEXP fitted_r, SEXP verbose_r, SEXP nReport_r);

  SEXP phenoMod0Pred(SEXP t_r, SEXP n_r, SEXP family_r, SEXP tNormBounds_r,
		     SEXP alpha_r, SEXP phi_r,
		     SEXP nSamples_r, SEXP fitted_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

}
