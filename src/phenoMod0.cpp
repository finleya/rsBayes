#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP phenoMod0(SEXP y_r, SEXP t_r, SEXP n_r, SEXP family_r, SEXP tNormBounds_r,
		 SEXP gamma_r,
                 SEXP phiIG_r, SEXP alphaUa_r, SEXP alphaUb_r,
                 SEXP phiStarting_r, SEXP alphaStarting_r,
                 SEXP phiTuning_r, SEXP alphaTuning_r,
                 SEXP nSamples_r, SEXP fitted_r, SEXP verbose_r, SEXP nReport_r){

    int i, s, inc = 1, nProtect=0;

    //get args
    double *y = REAL(y_r);
    double *t = REAL(t_r);
    int n = INTEGER(n_r)[0];

    int family = INTEGER(family_r)[0]; //0 is beta, 1 is normal, and 2 is truncated normal

    double tNormL, tNormU;
    if(family == 2){
      tNormL = REAL(tNormBounds_r)[0];
      tNormU = REAL(tNormBounds_r)[1];
    }

    double *gamma = REAL(gamma_r);
    
    double *phiIG = REAL(phiIG_r);
    double *alphaUa = REAL(alphaUa_r);
    double *alphaUb = REAL(alphaUb_r);

    double phiStarting = REAL(phiStarting_r)[0];
    double *alphaStarting = REAL(alphaStarting_r);

    double phiTuning = REAL(phiTuning_r)[0];
    double *alphaTuning = REAL(alphaTuning_r);

    int nSamples = INTEGER(nSamples_r)[0];
    int fitted = INTEGER(fitted_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
   
    int nAlpha = 7;

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      if(family == 0){
	Rprintf("Beta likelihood\n\n");
      }else if(family == 1){
	Rprintf("Normal likelihood\n\n");
      }else{
	Rprintf("Truncated Normal likelihood with bounds (%.3f, %.3f)\n\n", tNormL, tNormU);
      }
      
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Parameter,\tstarting,\ttuning,\t\tprior:\n");
      for(i = 0; i < nAlpha; i++){
	Rprintf("alpha.%i,\t%.5f,\t%.5f,\tUnif(%.5f, %.5f)\n", i+1, alphaStarting[i], alphaTuning[i], alphaUa[i], alphaUb[i]);
      }
      Rprintf("phi,\t\t%.5f,\t%.5f,\tIG(%.5f, %.5f)\n", i, phiStarting, phiTuning, phiIG[0], phiIG[1]);
    }

    //parameters
    int nTheta = nAlpha+1;//alpha, phi
    int alphaIndx = 0;
    int phiIndx = alphaIndx+nAlpha;

    //starting
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    for(i = 0; i < nTheta; i++){
      theta[i] = alphaStarting[i];
    }
    theta[phiIndx] = phiStarting;

    //tuning
    double *thetaTuning = (double *) R_alloc(nTheta, sizeof(double));
    for(i = 0; i < nAlpha; i++){
      thetaTuning[i] = alphaTuning[i];
    }
    thetaTuning[phiIndx] = phiTuning;

    //return stuff
    SEXP accept_r, thetaSamples_r, fittedSamples_r;
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++;

    if(fitted){
      PROTECT(fittedSamples_r = allocMatrix(REALSXP, n, nSamples)); nProtect++;
    }

    PROTECT(accept_r = allocVector(REALSXP, 2)); nProtect++;

    //other stuff
    double logPostCand, logPostCurrent;
    double *alpha = (double *) R_alloc(nAlpha, sizeof(double)); zeros(alpha, nAlpha);
    double *eta = (double *) R_alloc(n, sizeof(double));
    double phi;
    int accept = 0;
    int batchAccept = 0;
    int status = 0;
    double tmp;

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

      //current
      phi = theta[phiIndx];

      //likelihood
      G(theta, t, eta, n);

      if(family == 0){
	logPostCurrent = beta_logpost(n, y, eta, 1.0/phi);
      }else if(family == 1){
	logPostCurrent = norm_logpost(n, y, eta, phi);
      }else{
	logPostCurrent = tnorm_logpost(n, y, eta, phi, tNormL, tNormU);
      }

      logPostCurrent += -1.0*(1.0+phiIG[0])*log(phi)-phiIG[1]/phi+log(phi);

      //candidate
      alpha[2] = sampleU(theta[2], thetaTuning[2], alphaUa[2], alphaUb[2]);
      alpha[4] = sampleU(theta[4], thetaTuning[4], alphaUa[4], alphaUb[4]);
      alpha[5] = sampleU(theta[5], thetaTuning[5], alphaUa[5], alphaUb[5]);
      alpha[0] = sampleU(theta[0], thetaTuning[0], alphaUa[0], alphaUb[0]);
      alpha[1] = sampleU(theta[1], thetaTuning[1], alphaUa[1], gamma[1]-alpha[0]);      
      alpha[6] = sampleU(theta[6], thetaTuning[6], alphaUa[6], alphaUb[6]);
      alpha[3] = sampleU(theta[3], thetaTuning[3], alphaUa[3], alpha[6]);

      phi = exp(rnorm(log(theta[phiIndx]), thetaTuning[phiIndx]));

      //likelihood
      G(alpha, t, eta, n);

      if(family == 0){
	logPostCand = beta_logpost(n, y, eta, 1.0/phi);
      }else if(family == 1){
	logPostCand = norm_logpost(n, y, eta, phi);
      }else{
	logPostCand = tnorm_logpost(n, y, eta, phi, tNormL, tNormU);
      }

      logPostCand += -1.0*(1.0+phiIG[0])*log(phi)-phiIG[1]/phi+log(phi);

      //Metropolis accept/reject
      if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){
        F77_NAME(dcopy)(&nAlpha, alpha, &inc, theta, &inc);
        theta[phiIndx] = phi;
        accept++;
        batchAccept++;
      }

      ///////////////////////
      //save, fitted, report
      ///////////////////////
      F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[s*nTheta], &inc);

      if(fitted){
        G(theta, t, eta, n);
        phi = theta[phiIndx];

        for(i = 0; i < n; i++){
	  REAL(fittedSamples_r)[s*n+i] = eta[i];
        }
      }

      //report
      REAL(accept_r)[0] = 100.0*accept/nSamples;
      REAL(accept_r)[1] = 100.0*batchAccept/nReport;

      if(status == nReport){
        if(verbose){
      	Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
      	Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
      	Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/s);
      	Rprintf("-------------------------------------------------\n");
          #ifdef Win32
      	R_FlushConsole();
          #endif
        }
        batchAccept = 0;
        status = 0;
      }

      status++;

      R_CheckUserInterrupt();

    }//end sample loop

    PutRNGstate();

   //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 2;

    if(fitted){
      nResultListObjs++;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, thetaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.theta.samples"));

    SET_VECTOR_ELT(result_r, 1, accept_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("MH.acceptance"));

    if(fitted){
      SET_VECTOR_ELT(result_r, 2, fittedSamples_r);
      SET_VECTOR_ELT(resultName_r, 2, mkChar("p.fitted.samples"));
    }

    namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);

    return(result_r);
  }
}

