#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_bspline.h>
#include <omp.h>
#include <stdlib.h>

#include "../../libs/hbmTools.h"
#include "../../libs/kvpar.h"

#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <Rinternals.h>

using namespace std;

int main(int argc, char **argv) {
  
  string parfile;
  if (argc > 1) {
    parfile = argv[1];
  } else {
    parfile = "pfile-1";
  }
  
  /****************************************************************************
  * Obtain parameters from pfile
  ***************************************************************************/
  kvpar par(parfile);
  
  int seed; par.getVal("seed", seed);
  int nIter; par.getVal("n.iter", nIter);
  int nBatch; par.getVal("n.batch", nBatch);
  int nSamples = nIter * nBatch;
  string nChain; par.getVal("n.chain", nChain);
  int nThreads; par.getVal("n.threads", nThreads);
  string outFile; par.getVal("outfile", outFile);
  
  int nLocations; par.getVal("n.locations", nLocations);
  int nTime; par.getVal("n.time", nTime);
  int nDay; par.getVal("n.day", nDay);
  int nSpline; par.getVal("n.spline", nSpline);
  int nPhi; par.getVal("n.phi", nPhi);
  
  double *x; int n, a; par.getFile("x.file", x, n, a);
  if (n != (nTime * nDay * nLocations)) {
    cout << "DIMENSION PROBLEM" << endl;
  }
  double *y; int colY; par.getFile("y.file", y, n, colY);
  
  vector<double> betaTuning; par.getVal("beta.tuning", betaTuning);
  vector<double> betaStarting; par.getVal("beta.starting", betaStarting);
  vector<double> betaMu; par.getVal("beta.mu", betaMu);
  vector<double> betaSd; par.getVal("beta.sd", betaSd);
  
  double phiTuning; par.getVal("phi.tuning", phiTuning);
  double phiStarting; par.getVal("phi.starting", phiStarting);
  double phiA; par.getVal("phi.a", phiA);
  double phiB; par.getVal("phi.b", phiB);
  
  vector<double> lambdaStarting; par.getVal("lambda.starting", lambdaStarting);
  vector<double> lambdaTuning; par.getVal("lambda.tuning", lambdaTuning);
  double lambdaWiDf; par.getVal("lambda.df", lambdaWiDf);
  vector<double> lambdaIWSdiag; par.getVal("lambda.scale.diag", lambdaIWSdiag);
  
  double rhoTuning; par.getVal("rho.tuning", rhoTuning);
  double rhoStarting; par.getVal("rho.starting", rhoStarting);
  double rhoA; par.getVal("rho.a", rhoA);
  double rhoB; par.getVal("rho.b", rhoB);
  
  double wTune; par.getVal("w.tuning", wTune);
  
  set_seed(123, seed);
  omp_set_num_threads(nThreads);
  
  /****************************************************************************
  * Variables used in MCMC sampler
  ***************************************************************************/
  int i, j, b, c, d;
  char lower, transN, rside;
  int info, inc, nTimeLocDaySq, nTimeLocSq;
  int nTimeSq, nBreak, nLTR, nDaySq, nTimeLocDay, nTimeLoc, nLocSq;
  int nTheta, betaIndx, lambdaIndx, rhoIndx, phiIndx; 
  double roadMax, roadMin, mu, thetaCurr, logPostWCurr;
  double logPostThetaCurr, logPostThetaCand, rho, zero, one;
  double logPostWCand, wCurr, acceptRate, phi, *etaSamples;
  double *Z, *A, *w, *tmp, *wTuning, *eta, *wAccept, *lambda, *tmp_mm;
  double *theta, *thetaTuning, *accept, *B, *C, *lambdaIWS, *muSamples;
  double *thetaSamples, *wSamples, *fittedSamples, *zSamples, *Il;
  double detB, detL, detC, SLtrace, precision;
  
  /****************************************************************************
  * Log tuning parameters for use in adaptive MCMC sampler
  ***************************************************************************/
  nLTR = nDay * (nDay + 1) / 2;
  nTheta = nSpline + nPhi + nLTR + 1;
  betaIndx = 0;
  thetaTuning = new double[nTheta];
  phiIndx = betaIndx + nSpline;
  lambdaIndx = phiIndx + nPhi;
  rhoIndx = lambdaIndx + nLTR;
  for (i =0; i < nSpline; ++i) {
    thetaTuning[betaIndx + i] = log(betaTuning[i]);
  }
  thetaTuning[phiIndx] = log(phiTuning);
  for (i = 0; i < nLTR; ++i) {
    thetaTuning[lambdaIndx + i] = log(lambdaTuning[i]);
  }
  thetaTuning[rhoIndx] = log(rhoTuning);
  wTuning = new double[n];
  for (i = 0; i < n; ++i) {
    wTuning[i] = log(wTune);
  }
  
  /****************************************************************************
  * MCMC sampler preparation
  ***************************************************************************/
  theta = new double[nTheta];
  for (i = 0; i < nSpline; ++i) {
    theta[betaIndx + i] = betaStarting[i];
  }
  theta[phiIndx] = log(phiStarting);
  covTrans(lambdaStarting, &theta[lambdaIndx], nDay);
  theta[rhoIndx] = logit(rhoStarting, rhoA, rhoB);
  
  nTimeSq = nTime * nTime;
  nTimeLoc = nTime * nLocations;
  nLocSq = nLocations * nLocations;
  nTimeLocDay = nTime * nLocations * nDay;
  nTimeLocSq = nTimeLoc * nTimeLoc;
  nTimeLocDaySq = nTimeLocDay * nTimeLocDay;
  nDaySq = nDay * nDay;
  A = new double[nTimeSq];
  B = new double[nTimeLocSq];
  C = new double[nTimeLocDaySq];
  lower = 'L';
  info = 0;
  rside = 'R';
  w = new double[n]; zeros(w, n);
  lambda = new double[nDaySq];
  tmp_mm = new double[nDaySq];
  inc = 1;
  tmp = new double[n];
  zero = 0.0;
  one = 1.0;
  transN = 'N';
  accept = new double[nTheta]; zeros(accept, nTheta);
  eta = new double[n];
  wAccept = new double[n]; zeros(wAccept, n);
  thetaSamples = new double[nSamples * nTheta];
  wSamples = new double[nSamples * n];
  fittedSamples = new double[nSamples * n];
  muSamples = new double[nSamples * n];
  etaSamples = new double[nSamples * n];
  acceptRate = 0.43;
  Il = new double[nLocSq]; zeros(Il, nLocSq);
  for (i = 0; i < nLocations; ++i) {
    Il[i*nLocations + i] = 1.0;
  }
  lambdaIWS = new double[nDaySq]; zeros(lambdaIWS, nDaySq);
  for (i = 0; i < nDay; ++i) {
    lambdaIWS[i*nDay + i] = lambdaIWSdiag[i];
  }
  
  /****************************************************************************
  * Spline Coefficients
  ***************************************************************************/
  gsl_bspline_workspace *bw; 
  nBreak = nSpline - 2;
  bw = gsl_bspline_alloc(4, nBreak);
  gsl_vector *bc = gsl_vector_alloc(nSpline);
  roadMin = 1.0;
  roadMax = 0.0;
  for (i = 0; i < n * colY; ++i) {
    if (x[i] < roadMin) {
      roadMin = x[i];
    }
    if (x[i] > roadMax) {
      roadMax = x[i];
    }
  }
  // gsl_bspline_knots_uniform(roadMin, roadMax, bw);
  gsl_vector *knots = gsl_vector_alloc(nBreak);
  gsl_vector_set(knots, 0, roadMin);
  gsl_vector_set(knots, 1, 1);
  gsl_vector_set(knots, 2, roadMax);
  gsl_bspline_knots(knots, bw);
  Z = new double[n * nSpline];
  zSamples = new double[nSamples * n * nSpline];
  
  for (i = 0; i < n; ++i) {
    gsl_bspline_eval(x[i], bc, bw);
    for (j = 0; j < nSpline; ++j) {
      Z[j * n + i] = gsl_vector_get(bc, j);
    }
  }
  
  /****************************************************************************
  * Sampler
  ***************************************************************************/
  
  for (a = 0, d = 0; a < nIter; ++a) {
    for (b = 0; b < nBatch; ++b, ++d) {

      /************************************
      * Current theta
      ***********************************/
      logPostThetaCurr = 0.0;
      phi = exp(theta[phiIndx]);
      covTransInvExpand(&theta[lambdaIndx], lambda, nDay);
      rho = logitInv(theta[rhoIndx], rhoA, rhoB);
      
      /************
      * Priors
      ***********/
      for (i = 0; i < nSpline; ++i) {
        logPostThetaCurr += dnorm(theta[betaIndx + i], betaMu[i], betaSd[i], 1);
      }
      logPostThetaCurr += iGammaLogpostJacobian(phi, phiA, phiB);
      logPostThetaCurr += uniformLogpostJacobian(rho, rhoA, rhoB);
      
      detL = 0.0;
      SLtrace = 0.0;
      for (i = 0; i < nDay; ++i) {
        detL += 2 * log(lambda[i*nDay + i]);
      }
      for(i = 0; i < nDay; ++i){
        logPostThetaCurr += (nDay-i)*log(lambda[i*nDay+i])+log(lambda[i*nDay+i]);
      }
      //get S*L^-1, already have the chol of L (i.e., lambda)
      dpotri_(&lower, &nDay, lambda, &nDay, &info); 
      if(info != 0) {
        cout << "c++ error: dpotrf failed" << endl;
      }
      dsymm_(&rside, &lower, &nDay, &nDay, &one, lambda, &nDay, lambdaIWS, &nDay, &zero, tmp_mm, &nDay);
      for(i = 0; i < nDay; ++i){
        SLtrace += tmp_mm[i*nDay+i];
      }
      logPostThetaCurr += -0.5*(lambdaWiDf+nDay+1)*detL - 0.5*SLtrace;

      AR1(nTime, rho, one, A);
      detB = 0.0;
      dpotrf_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotrf failed" << endl;
      }
      for (i = 0; i < nTime; ++i) {
        detB += 2 * log(A[i*nTime + i]);
      }
      detB *= nLocations;
      dpotri_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotri failed" << endl;
      }
      fillUTri(A, nTime);
      kron(Il, nLocations, nLocations, A, nTime, nTime, B, nTimeLoc, nTimeLoc);
      fillUTri(lambda, nDay);
      kron(lambda, nDay, nDay, B, nTimeLoc, nTimeLoc, C, nTimeLocDay, nTimeLocDay);
      detC = nTimeLoc * detL + nDay * detB;
      dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
      logPostThetaCurr += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
      
      /************
      * Likelihood
      ***********/
      dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
             eta, &inc);
      daxpy_(&n, &one, w, &inc, eta, &inc);
      logPostThetaCurr += betaLogpostJacobian(n, y, eta, phi);
      
      /************************************
      * Proposed theta
      ***********************************/
      for (c = 0; c < nTheta; ++c) {
        thetaCurr = theta[c];
        theta[c] = rnorm(thetaCurr, exp(thetaTuning[c]));
        logPostThetaCand = 0.0;
        phi = exp(theta[phiIndx]);
        covTransInvExpand(&theta[lambdaIndx], lambda, nDay);
        rho = logitInv(theta[rhoIndx], rhoA, rhoB);
        
        /************
        * Priors
        ***********/
        for (i = 0; i < nSpline; ++i) {
          logPostThetaCand += dnorm(theta[betaIndx + i], betaMu[i], betaSd[i], 1);
        }
        logPostThetaCand += iGammaLogpostJacobian(phi, phiA, phiB);
        logPostThetaCand += uniformLogpostJacobian(rho, rhoA, rhoB);
        
        detL = 0.0;
        SLtrace = 0.0;
        for (i = 0; i < nDay; ++i) {
          detL += 2 * log(lambda[i*nDay + i]);
        }
        for(i = 0; i < nDay; ++i){
          logPostThetaCand += (nDay-i)*log(lambda[i*nDay+i])+log(lambda[i*nDay+i]);
        }
        //get S*L^-1, already have the chol of L (i.e., lambda)
        dpotri_(&lower, &nDay, lambda, &nDay, &info); 
        if(info != 0) {
          cout << "c++ error: dpotrf failed" << endl;
        }
        dsymm_(&rside, &lower, &nDay, &nDay, &one, lambda, &nDay, lambdaIWS, &nDay, &zero, tmp_mm, &nDay);
        for(i = 0; i < nDay; ++i){
          SLtrace += tmp_mm[i*nDay+i];
        }
        logPostThetaCand += -0.5*(lambdaWiDf+nDay+1)*detL - 0.5*SLtrace;
        
        AR1(nTime, rho, one, A);
        detB = 0.0;
        dpotrf_(&lower, &nTime, A, &nTime, &info);
        if (info) {
          cout << "dpotrf failed" << endl;
        }
        for (i = 0; i < nTime; ++i) {
          detB += 2 * log(A[i*nTime + i]);
        }
        detB *= nLocations;
        dpotri_(&lower, &nTime, A, &nTime, &info);
        if (info) {
          cout << "dpotri failed" << endl;
        }
        fillUTri(A, nTime);
        kron(Il, nLocations, nLocations, A, nTime, nTime, B, nTimeLoc, nTimeLoc);
        fillUTri(lambda, nDay);
        kron(lambda, nDay, nDay, B, nTimeLoc, nTimeLoc, C, nTimeLocDay, nTimeLocDay);
        detC = nTimeLoc * detL + nDay * detB;
        dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
        logPostThetaCand += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
        
        /************
        * Likelihood
        ***********/
        dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
               eta, &inc);
        daxpy_(&n, &one, w, &inc, eta, &inc);
        logPostThetaCand += betaLogpostJacobian(n, y, eta, phi);
        
        /************
        * Accept/Reject
        ***********/
        if (runif(0.0, 1.0) <= exp(logPostThetaCand - logPostThetaCurr)) {
	        logPostThetaCurr = logPostThetaCand;
	        accept[c]++;
        } else {
          theta[c] = thetaCurr;
        }
      }
      
      /************************************
      * Current random effects
      ***********************************/
      
      logPostWCurr = 0.0;
      phi = exp(theta[phiIndx]);
      covTransInvExpand(&theta[lambdaIndx], lambda, nDay);
      rho = logitInv(theta[rhoIndx], rhoA, rhoB);
      
      detL = 0.0;
      for (i = 0; i < nDay; ++i) {
        detL += 2 * log(lambda[i*nDay + i]);
      }
      dpotri_(&lower, &nDay, lambda, &nDay, &info); 
      if(info != 0) {
        cout << "c++ error: dpotrf failed" << endl;
      }
      AR1(nTime, rho, one, A);
      detB = 0.0;
      dpotrf_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotrf failed" << endl;
      }
      for (i = 0; i < nTime; ++i) {
        detB += 2 * log(A[i*nTime + i]);
      }
      detB *= nLocations;
      dpotri_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotri failed" << endl;
      }
      fillUTri(A, nTime);
      kron(Il, nLocations, nLocations, A, nTime, nTime, B, nTimeLoc, nTimeLoc);
      fillUTri(lambda, nDay);
      kron(lambda, nDay, nDay, B, nTimeLoc, nTimeLoc, C, nTimeLocDay, nTimeLocDay);
      detC = nTimeLoc * detL + nDay * detB;
      dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
      logPostWCurr += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
      
      /************
      * Likelihood
      ***********/
      dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
             eta, &inc);
      dcopy_(&n, eta, &inc, tmp, &inc);
      daxpy_(&n, &one, w, &inc, tmp, &inc);
      logPostWCurr += betaLogpostJacobian(n, y, tmp, phi);
      
      /************************************
      * Proposed random effects
      ***********************************/
      for (c = 0; c < n; ++c) {
        wCurr = w[c];
        w[c] = rnorm(wCurr, exp(wTuning[c]));
        
        logPostWCand = 0.0;
        
        dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
        logPostWCand += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
        
        /************
        * Likelihood
        ***********/
        dcopy_(&n, eta, &inc, tmp, &inc);
        daxpy_(&n, &one, w, &inc, tmp, &inc);
        logPostWCand += betaLogpostJacobian(n, y, tmp, phi);
        
        if (runif(0.0, 1.0) <= exp(logPostWCand - logPostWCurr)) {
          logPostWCurr = logPostWCand;
	        wAccept[c]++;
        } else {
          w[c] = wCurr;
        }
      }
      
      /************************************
      * Save current values
      ***********************************/
      dcopy_(&nTheta, theta, &inc, &thetaSamples[nTheta * d], &inc);
      dcopy_(&n, w ,&inc, &wSamples[n * d], &inc);
      
      dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
             eta, &inc);

      daxpy_(&n, &one, w, &inc, eta, &inc);
      
      phi = exp(theta[phiIndx]);

      for (i = 0; i < n; ++i) {
        mu = logitInv(eta[i], 0.0, 1.0);
        muSamples[d*n + i] = mu;
        fittedSamples[d * n + i] = rbeta(mu * phi, (1.0 - mu) * phi);
      }
    } // batch
    
    cout << "-------------------------------" << endl;
    
    cout << "percent complete: " << 100 * d / nSamples << endl;
    
    /************************************
    * Adjust tuning parameters
    ***********************************/
    for (i = 0; i < nTheta; ++i) {
      cout << "theta[" << i << "] acceptance: " << 100 * accept[i] / nBatch <<
        " tuning: " << exp(thetaTuning[i]) << endl;
      if (static_cast<double> (accept[i]) / nBatch > acceptRate) {
        thetaTuning[i] += min(0.01, 1.0 / sqrt(static_cast<double> (a)));
      } else {
        thetaTuning[i] -= min(0.01, 1.0 / sqrt(static_cast<double> (a)));
      }
      accept[i] = 0;
    }
    
    for (i = 0; i < n; ++i) {
      if (i < 5) {
        cout << "w[" << i << "] acceptance: " << 100 * wAccept[i] / nBatch <<
          " tuning: " << exp(wTuning[i]) << endl;
      }
      if (static_cast<double> (wAccept[i]) / nBatch > acceptRate) {
        wTuning[i] += min(0.01, 1.0 / sqrt(static_cast<double> (a)));
      } else {
        wTuning[i] -= min(0.01, 1.0 / sqrt(static_cast<double> (a)));
      }
      wAccept[i] = 0;
    }
  } // iterations
  
  /****************************************************************************
   * Write output
   ***************************************************************************/
  for (i = 0; i < nSamples; ++i) {
    thetaSamples[i*nTheta + phiIndx] = exp(thetaSamples[i*nTheta + phiIndx]);
    covTransInv(&thetaSamples[nTheta*i + lambdaIndx], &thetaSamples[nTheta*i + lambdaIndx], 
                nDay);
    thetaSamples[i*nTheta + rhoIndx] = logitInv(thetaSamples[i*nTheta + rhoIndx],
                                                rhoA, rhoB);
  }
  
  writeRMatrix(outFile + "-fitted-" + nChain, fittedSamples, n, nSamples);
  writeRMatrix(outFile + "-muSamples-" + nChain, muSamples, n, nSamples);
  writeRMatrix(outFile + "-theta-" + nChain, thetaSamples, nTheta, nSamples);
  writeRMatrix(outFile + "-w-" + nChain, wSamples, n, nSamples);
  writeRMatrix(outFile + "-spline-setup", Z, n, nSpline);
  
  return 0;
}
