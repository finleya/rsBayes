#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif


extern "C" {
  
  SEXP phenoMod0Pred(SEXP t_r, SEXP n_r, SEXP family_r, SEXP tNormBounds_r,
		     SEXP alpha_r, SEXP phi_r,
		     SEXP nSamples_r, SEXP fitted_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){

    int i, s, inc = 1, nProtect=0;
    
    //get args
    double *t = REAL(t_r);
    int n = INTEGER(n_r)[0];
    int family = INTEGER(family_r)[0]; //0 is beta, 1 is normal, and 2 is truncated normal
    double tNormL, tNormU;
    if(family == 2){
      tNormL = REAL(tNormBounds_r)[0];
      tNormU = REAL(tNormBounds_r)[1];
    }
    
    double *alpha = REAL(alpha_r);
    double *phi = REAL(phi_r);
    int nSamples = INTEGER(nSamples_r)[0];
    int fitted = INTEGER(fitted_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
    int nAlpha = 7;
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      if(fitted == 1){
	Rprintf("\tSampling fitted values\n");
      }else{
	Rprintf("\tSampling predicted values\n");
      }
      Rprintf("----------------------------------------\n");
      Rprintf("Sampling at %i time points\n", n);
      if(family == 0){
	Rprintf("Beta likelihood\n\n");
      }else if(family == 1){
	Rprintf("Normal likelihood\n\n");
      }else{
	Rprintf("Truncated Normal likelihood with bounds (%.3f, %.3f)\n\n", tNormL, tNormU);
      }
      
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      #ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    }

    //return stuff
    SEXP ySamples_r;
    PROTECT(ySamples_r = allocMatrix(REALSXP, n, nSamples)); nProtect++;
    
    //other stuff
    double logPostCand, logPostCurrent;
    int accept = 0;
    int batchAccept = 0;
    int status = 0;
     
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("----------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    GetRNGstate();

    for(s = 0; s < nSamples; s++){
     
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(i = 0; i < n; i++){
	REAL(ySamples_r)[s*n+i] = G(&alpha[s*nAlpha], t[i]);
      }

      //add prediction variance
      if(fitted == 0){
	for(i = 0; i < n; i++){
	  if(family == 0){
	    REAL(ySamples_r)[s*n+i] = rbeta(REAL(ySamples_r)[s*n+i]*1.0/phi[s], (1.0-REAL(ySamples_r)[s*n+i])*1.0/phi[s]);
	  }else if(family == 1){
	    REAL(ySamples_r)[s*n+i] = rnorm(REAL(ySamples_r)[s*n+i], sqrt(phi[s]));
	  }else{
	    REAL(ySamples_r)[s*n+i] = rtnorm(REAL(ySamples_r)[s*n+i], sqrt(phi[s]), tNormL, tNormU);
	  }
	}
      }

      if(status == nReport){
        if(verbose){
      	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
      	  Rprintf("-------------------------------------------------\n");
#ifdef Win32
      	  R_FlushConsole();
#endif
        }
        status = 0;
      }
      
      status++;
      
      R_CheckUserInterrupt();

    }//end sample loop
    
    PutRNGstate();

   //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 1;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, ySamples_r);
    if(fitted == 0){
      SET_VECTOR_ELT(resultName_r, 0, mkChar("p.predictive.samples"));
    }else{
      SET_VECTOR_ELT(resultName_r, 0, mkChar("p.fitted.samples"));
    }
    
    namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);

    return(result_r);
  }
}

