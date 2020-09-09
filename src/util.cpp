#include <string>
#include <limits>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "dtnorm_data.h"

double beta_logpost(int &n, double *Y, double *eta, double phi){

  int i;
  double loglike = 0.0;

  for(i = 0; i < n; i++){
    loglike += dbeta(Y[i], eta[i]*phi, (1.0-eta[i])*phi, true);   
    //loglike += lgammafn(phi) - lgammafn(mu*phi) - lgammafn((1.0-mu)*phi) + (mu*phi - 1.0)*log(Y[i]) + ((1.0 - mu)*phi - 1.0)*log(1.0 - Y[i]);
    //loglike += log(gamma(phi)) - log(gamma(mu*phi)) - log(gamma((1.0-mu)*phi)) + (mu*phi - 1.0)*log(Y[i]) + ((1.0 - mu)*phi - 1.0)*log(1.0 - Y[i]);
  }

  return loglike;  
}

double norm_logpost(int &n, double *Y, double *eta, double phi){

  int i;
  double loglike = 0.0;

  for(i = 0; i < n; i++){
    loglike += dnorm(eta[i], Y[i], sqrt(phi), true);
  }

  return loglike;  
}


double tnorm_logpost(int &n, double *Y, double *eta, double phi, double a, double b){

  int i;
  double loglike = 0.0;

  for(i = 0; i < n; i++){
    loglike += dtnorm(eta[i], Y[i], sqrt(phi), a, b, 1);
  }

  return loglike;  
}



void G(double *alpha, double *t, double *eta, int n){
  
  int i;
  double delta = (alpha[2]*alpha[3]+alpha[5]*alpha[6])/(alpha[2]+alpha[5]);
  
  for(i = 0; i < n; i++){
    
    if(t[i] <= delta){
      
      eta[i] = alpha[0] + (alpha[1] - alpha[4]*t[i])/(1.0 + exp(-alpha[2]*(t[i]-alpha[3])));
      
    }else{
      
      eta[i] = alpha[0] + (alpha[1] - alpha[4]*t[i])/(1.0 + exp(-alpha[5]*(alpha[6]-t[i])));
      
    }
    
  }
  
}

double G(double *alpha, double t){
  
  if(t <= (alpha[2]*alpha[3]+alpha[5]*alpha[6])/(alpha[2]+alpha[5])){
    
    return alpha[0] + (alpha[1] - alpha[4]*t)/(1.0 + exp(-alpha[2]*(t-alpha[3])));
    
  }else{
    
    return alpha[0] + (alpha[1] - alpha[4]*t)/(1.0 + exp(-alpha[5]*(alpha[6]-t)));
    
  }
  
}

double mirror(double theta, double a, double b){
 
  if(theta > b){     
    return mirror(theta - 2.0*fabs(theta-b), a, b);
  }else if(theta < a){
    return mirror(theta + 2.0*fabs(a-theta), a, b);  
  }else{
    return theta;
  }
}

double sampleU(double mu, double tuning, double a, double b){
  return mirror(rnorm(mu, tuning), a, b);
}


void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}


void printMtrx(double *m, int nRow, int nCol){

  int i, j;

  for(i = 0; i < nRow; i++){
    Rprintf("\t");
    for(j = 0; j < nCol; j++){
      Rprintf("%.10f\t", m[j*nRow+i]);
    }
    Rprintf("\n");
  }
}


double dtnorm(double x, double mu, double sigma, double a, double b, int lg){
    
  double denom, xtmp;
  
  if((x < a) | (x > b)){
    return(R_NegInf);
  }
  
  if((x >= a) & (x <= b)){  
    denom = pnorm(b, mu, sigma, 1, 0) - pnorm(a, mu, sigma, 1, 0);
    xtmp = dnorm(x, mu, sigma, lg);
    
    if(lg){
      xtmp = xtmp - log(denom);
    }else{
      xtmp = xtmp/denom;
    }
  }
  
  return(xtmp);  
}

// Andrew Finley modified Alan R. Rogers's code to use R's rngs. Robert's code
// is at https://github.com/alanrogers/dtnorm Rogers code was based on C++
// original written by G. Dolle and V Mazet, which is available at
// http://miv.u-strasbg.fr/mazet/rtnorm/rtnormCpp.zip.  That original code
// implemented several methods, including Chopin's algorithm (2011. Stat Comput. 21:275-288).
// Different methods were used in different regions of parameter space. 

// Roberts notes own numerical experiments indicated that this method can be
// improved by changing the criteria used in selecting algorithms. In
// some regions of parameter space, it is best to use Robert's method
// (1995. Statistics and Computing, 5(2):121). The current code
// implements the revised method in C.

// Licence: GNU General Public License Version 2
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version. This program is distributed in the hope that
// it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details. You should have received a
// copy of the GNU General Public License along with this program; if not,
// see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt

double yl(int k) {
  int N = 4001;           // Index of the right tail
  double      yl0 = 0.053513975472;   // y_l of the leftmost rectangle
  double      ylN = 0.000914116389555;    // y_l of the rightmost rectangle
  
  if(k == 0)
    return yl0;
  
  else if(k == N - 1)
    return ylN;
  
  else if(k <= 1953)
    return yu[k - 1];
  
  else
    return yu[k + 1];
}

/// Rejection algorithm with a truncated exponential proposal
double rtexp(double a, double b) {
    double      twoasq = 2*a*a;
    double      expab = expm1(-a * (b - a));
    double      z, e;

    do{
      //z = log1p(gsl_rng_uniform(rng) * expab);
      //e = -log(gsl_rng_uniform(rng));
      z = log1p(runif(0.0, 1.0) * expab);
      e = -log(runif(0.0, 1.0));
      
    }while(twoasq*e <= z*z);
    return a - z/a;
}

double rtnorm(const double mu, const double sigma, double a, double b) {
    // Design variables
  int N = 4001;           // Index of the right tail
  //double      xmin = -2.00443204036;  // Left bound
  double      xmax = 3.48672170399;   // Right bound
  int         kmin = 5;
  double      INVH = 1631.73284006;   // = 1/h, h being the minimal interval range
  int         I0 = 3271;      // = - floor(x(0)/h)
  double      ALPHA = 1.837877066409345;  // = log(2*pi)
  int         xsize = sizeof(x) / sizeof(double); // Length of table x
  int         stop = false;  
  double      r, z, e, ylk, simy, lbound, u, d, sim;
  int         i, ka, kb, k;
  
  // Scaling
  if(mu != 0 || sigma != 1) {
    a = (a - mu) / sigma;
    b = (b - mu) / sigma;
  }
  
  double bma = b - a;
  
  // Check if a < b
  if(a >= b) {
    error("bound b must be greater than bound a!\n");
  }
  
  // Check if |a| < |b|
  else if(fabs(a) > fabs(b))
    r = -rtnorm(0.0, 1.0, -b, -a);
  
  // Truncated exponential proposal
  else if(a > 3.0 && bma > 0.5)
    r = rtexp(a, b);
  
  // Gaussian proposal
  else if(a < 0.0 && b > 0.0 && bma > 1.0) {
    while(!stop) {
      r = rnorm(0, 1.0);
      stop = (r >= a) && (r <= b);
    }
  }
  
  // Robert's algorithm
  else if(b < 0.0
	  || (a < 0.0 && (b < 0.0 || bma < 1.0))
	  || (a > 3.0 && bma < 0.5)) {
    double rho;
    do{
      r = runif(a, b);
      if( a*b < 0.0 )               // of opposite sign
	rho = exp(-0.5*r*r);
      else if( b < 0.0 )             // both negative
	rho = exp( 0.5*(b*b - r*r) );
      else{                           // both positive
	//assert(a > 0.0);
	if(!(a > 0.0)){error("error in Robert's algorithm\n");}
	rho = exp( 0.5*(a*a - r*r) );
      }
      
      u = runif(0.0, 1.0);
    }while( u > rho );
  }
  
  // Chopin's algorithm
  else {
    if(!(0.0 <= a && a <= 3.0)) {
      error("error in Chopin's algorithm\n");
    }
    if(!(0.0 <= a)){error("error in Chopin's algorithm\n");}
    if(!(a <= 3.0)){error("error in Chopin's algorithm\n");}
    
    // Compute ka
    i = I0 + floor(a * INVH);
    ka = ncell[i];
    
    // Compute kb
    (b >= xmax) ?
      kb = N : (i = I0 + floor(b * INVH), kb = ncell[i]
		);
    
    // If |b-a| is small, use rejection algorithm with a truncated
    // exponential proposal
    if(abs(kb - ka) < kmin) {
      r = rtexp(a, b);
      stop = true;
    }
    
    while(!stop) {
      // Sample integer between ka and kb
      k = floor(runif(0.0, 1.0) * (kb - ka + 1)) + ka;
      
      if(k == N) {
	// Right tail
	lbound = x[xsize - 1];
	z = -log(runif(0.0, 1.0));
	e = -log(runif(0.0, 1.0));
	z = z / lbound;
	
	if((z*z <= 2 * e) && (z < b - lbound)) {
	  // Accept this proposition, otherwise reject
	  r = lbound + z;
	  stop = true;
	}
      }
      
      else if((k <= ka + 1) || (k >= kb - 1 && b < xmax)) {
	
	// Two leftmost and rightmost regions
	sim = x[k] + (x[k + 1] - x[k]) * runif(0.0, 1.0);
	
	if((sim >= a) && (sim <= b)) {
	  // Accept this proposition, otherwise reject
	  simy = yu[k] * runif(0.0, 1.0);
	  if((simy < yl(k))
	     || (sim * sim + 2 * log(simy) + ALPHA) < 0) {
	    r = sim;
	    stop = true;
	  }
	}
      }
      
      else                // All the other boxes
	{
	  u = runif(0.0, 1.0);
	  simy = yu[k] * u;
	  d = x[k + 1] - x[k];
	  ylk = yl(k);
	  if(simy < ylk)  // That's what happens most of the time
	    {
	      r = x[k] + u * d * yu[k] / ylk;
	      stop = true;
	    } else {
	    sim = x[k] + d * runif(0.0, 1.0);
	    
	    // Otherwise, check you're below the pdf curve
	    if((sim * sim + 2 * log(simy) + ALPHA) < 0) {
	      r = sim;
	      stop = true;
	    }
	  }
	  
	}
    }
  }
  
  // Scaling
  if(mu != 0 || sigma != 1)
    r = r * sigma + mu;
  
  return r;
}
