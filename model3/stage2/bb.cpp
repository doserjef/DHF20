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
  double *y; int colY; par.getFile("y.file", y, n, colY);
  
  vector<double> betaTuning; par.getVal("beta.tuning", betaTuning);
  vector<double> betaStarting; par.getVal("beta.starting", betaStarting);
  vector<double> betaMu; par.getVal("beta.mu", betaMu);
  vector<double> betaSd; par.getVal("beta.sd", betaSd);
  
  vector<double> phiStarting; par.getVal("phi.starting", phiStarting);
  vector<double> phiTuning; par.getVal("phi.tuning", phiTuning);
  vector<double> phiA; par.getVal("phi.a", phiA);
  vector<double> phiB; par.getVal("phi.b", phiB);
  
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
  double alphaMax, alphaMin, mu, thetaCurr, logPostWCurr;
  double logPostThetaCurr, logPostThetaCand, rho, zero, one;
  double logPostWCand, wCurr, acceptRate;
  double *Z, *A, *w, *tmp, *wTuning, *eta, *wAccept, *lambda, *tmp_mm;
  double *theta, *thetaTuning, *accept, *B, *C, *lambdaIWS;
  double *thetaSamples, *wSamples, *fittedSamples, *Il;
  double detB, detL, detC, SLtrace;
  double *phi, *varDesign, *sigma, precision;
  
  /****************************************************************************
  * Log tuning parameters for use in adaptive MCMC sampler
  ***************************************************************************/
  nLTR = nDay * (nDay + 1) / 2;
  // nTheta = nSpline + nPhi + nLTR + 1;
  nTheta = nSpline + nPhi;
  betaIndx = 0;
  thetaTuning = new double[nTheta];
  phiIndx = betaIndx + nSpline;
  lambdaIndx = phiIndx + nPhi;
  rhoIndx = lambdaIndx + nLTR;
  for (i =0; i < nSpline; ++i) {
    thetaTuning[betaIndx + i] = log(betaTuning[i]);
  }
  for (i = 0; i < nPhi; ++i) {
    thetaTuning[phiIndx + i] = log(phiTuning[i]);
  }
  // for (i = 0; i < nLTR; ++i) {
  //   thetaTuning[lambdaIndx + i] = log(lambdaTuning[i]);
  // }
  // thetaTuning[rhoIndx] = log(rhoTuning);
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
  for (i = 0; i < nPhi; ++i) {
    theta[phiIndx + i] = log(phiStarting[i]);
  }
  // covTrans(lambdaStarting, &theta[lambdaIndx], nDay);
  // theta[rhoIndx] = logit(rhoStarting, rhoA, rhoB);
  
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
  fittedSamples = new double[nSamples *n];
  acceptRate = 0.43;
  Il = new double[nLocSq]; zeros(Il, nLocSq);
  for (i = 0; i < nLocations; ++i) {
    Il[i*nLocations + i] = 1.0;
  }
  lambdaIWS = new double[nDaySq]; zeros(lambdaIWS, nDaySq);
  for (i = 0; i < nDay; ++i) {
    lambdaIWS[i*nDay + i] = lambdaIWSdiag[i];
  }
  phi = new double[nPhi];
  sigma = new double[n];
  varDesign = new double[n * nPhi]; ones(varDesign, n);
  
  /****************************************************************************
  * Spline Coefficients
  ***************************************************************************/
  gsl_bspline_workspace *bw; 
  nBreak = nSpline - 2;
  bw = gsl_bspline_alloc(4, nBreak);
  gsl_vector *bc = gsl_vector_alloc(nSpline);
  alphaMin = 1.0;
  alphaMax = 0.0;
  for (i = 0; i < n * colY; ++i) {
    if (x[i] < alphaMin) {
      alphaMin = x[i];
    }
    if (x[i] > alphaMax) {
      alphaMax = x[i];
    }
  }
  gsl_vector *knots = gsl_vector_alloc(nBreak);
  gsl_vector_set(knots, 0, 0.0);
  gsl_vector_set(knots, 1, 0.7);
  gsl_vector_set(knots, 2, 0.8);
  gsl_vector_set(knots, 3, 0.97);
  gsl_vector_set(knots, 4, 0.99);
  gsl_vector_set(knots, 5, 1.0);
  gsl_bspline_knots(knots, bw);
  Z = new double[n * nSpline];
  
  // for (i = 0; i < n; ++i) {
  //   gsl_bspline_eval(x[i], bc, bw);
  //   for (j = 0; j < nSpline; ++j) {
  //     Z[j * n + i] = gsl_vector_get(bc, j);
  //   }
  // }
  
  
  /****************************************************************************
  * Sampler
  ***************************************************************************/
  
  for (a = 0, d = 0; a < nIter; ++a) {
    for (b = 0; b < nBatch; ++b, ++d) {
      
      /************************************
       * Update spline matrix and save
       ***********************************/
      for (i = 0; i < n; ++i) {
        gsl_bspline_eval(x[d*n + i], bc, bw);
        for (j = 0; j < nSpline; ++j) {
          Z[j * n + i] = gsl_vector_get(bc, j);
        }
        varDesign[n + i] = exp(x[d*n + i]);
        // varDesign[n + i] = pow(x[d*n +i], 2);
      }
      
      
      /************************************
      * Current theta
      ***********************************/
      logPostThetaCurr = 0.0;
      for (i = 0; i < nPhi; ++i) {
        phi[i] = exp(theta[phiIndx + i]);
      }
      // covTransInvExpand(&theta[lambdaIndx], lambda, nDay);
      // rho = logitInv(theta[rhoIndx], rhoA, rhoB);
      
      /************
      * Priors
      ***********/
      for (i = 0; i < nSpline; ++i) {
        logPostThetaCurr += dnorm(theta[betaIndx + i], betaMu[i], betaSd[i], 1);
      }
      for (i = 0; i < nPhi; ++i) {
        logPostThetaCurr += uniformLogpostJacobian(phi[i], phiA[i], phiB[i]);
      }
      // logPostThetaCurr += uniformLogpostJacobian(rho, rhoA, rhoB);
      // 
      // detL = 0.0;
      // SLtrace = 0.0;
      // for (i = 0; i < nDay; ++i) {
      //   detL += 2 * log(lambda[i*nDay + i]);
      // }
      // for(i = 0; i < nDay; ++i){
      //   logPostThetaCurr += (nDay-i)*log(lambda[i*nDay+i])+log(lambda[i*nDay+i]);
      // }
      // //get S*L^-1, already have the chol of L (i.e., lambda)
      // dpotri_(&lower, &nDay, lambda, &nDay, &info);
      // if(info != 0) {
      //   cout << "c++ error: dpotrf failed" << endl;
      // }
      // dsymm_(&rside, &lower, &nDay, &nDay, &one, lambda, &nDay, lambdaIWS, &nDay, &zero, tmp_mm, &nDay);
      // for(i = 0; i < nDay; ++i){
      //   SLtrace += tmp_mm[i*nDay+i];
      // }
      // logPostThetaCurr += -0.5*(lambdaWiDf+nDay+1)*detL - 0.5*SLtrace;
      // 
      // AR1(nTime, rho, one, A);
      // detB = 0.0;
      // dpotrf_(&lower, &nTime, A, &nTime, &info);
      // if (info) {
      //   cout << "dpotrf failed" << endl;
      // }
      // for (i = 0; i < nTime; ++i) {
      //   detB += 2 * log(A[i*nTime + i]);
      // }
      // detB *= nLocations;
      // dpotri_(&lower, &nTime, A, &nTime, &info);
      // if (info) {
      //   cout << "dpotri failed" << endl;
      // }
      // fillUTri(A, nTime);
      // kron(Il, nLocations, nLocations, A, nTime, nTime, B, nTimeLoc, nTimeLoc);
      // fillUTri(lambda, nDay);
      // kron(lambda, nDay, nDay, B, nTimeLoc, nTimeLoc, C, nTimeLocDay, nTimeLocDay);
      // detC = nTimeLoc * detL + nDay * detB;
      // dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
      // logPostThetaCurr += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
      
      /************
      * Likelihood
      ***********/
      dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
             eta, &inc);
      // daxpy_(&n, &one, w, &inc, eta, &inc);
      dgemv_(&transN, &n, &nPhi, &one, varDesign, &n, phi, &inc, &zero, 
             sigma, &inc);
      logPostThetaCurr += betaLogpostJacobianVar(n, y, eta, sigma);
      
      /************************************
      * Proposed theta
      ***********************************/
      for (c = 0; c < nTheta; ++c) {
        thetaCurr = theta[c];
        theta[c] = rnorm(thetaCurr, exp(thetaTuning[c]));
        logPostThetaCand = 0.0;
        for (i = 0; i < nPhi; ++i) {
          phi[i] = exp(theta[phiIndx + i]);
        }
        // covTransInvExpand(&theta[lambdaIndx], lambda, nDay);
        // rho = logitInv(theta[rhoIndx], rhoA, rhoB);
        
        /************
        * Priors
        ***********/
        for (i = 0; i < nSpline; ++i) {
          logPostThetaCand += dnorm(theta[betaIndx + i], betaMu[i], betaSd[i], 1);
        }
        for (i = 0; i < nPhi; ++i) {
          logPostThetaCand += uniformLogpostJacobian(phi[i], phiA[i], phiB[i]);
        }
        // logPostThetaCand += uniformLogpostJacobian(rho, rhoA, rhoB);
        // 
        // detL = 0.0;
        // SLtrace = 0.0;
        // for (i = 0; i < nDay; ++i) {
        //   detL += 2 * log(lambda[i*nDay + i]);
        // }
        // for(i = 0; i < nDay; ++i){
        //   logPostThetaCand += (nDay-i)*log(lambda[i*nDay+i])+log(lambda[i*nDay+i]);
        // }
        // //get S*L^-1, already have the chol of L (i.e., lambda)
        // dpotri_(&lower, &nDay, lambda, &nDay, &info);
        // if(info != 0) {
        //   cout << "c++ error: dpotrf failed" << endl;
        // }
        // dsymm_(&rside, &lower, &nDay, &nDay, &one, lambda, &nDay, lambdaIWS, &nDay, &zero, tmp_mm, &nDay);
        // for(i = 0; i < nDay; ++i){
        //   SLtrace += tmp_mm[i*nDay+i];
        // }
        // logPostThetaCand += -0.5*(lambdaWiDf+nDay+1)*detL - 0.5*SLtrace;
        // 
        // AR1(nTime, rho, one, A);
        // detB = 0.0;
        // dpotrf_(&lower, &nTime, A, &nTime, &info);
        // if (info) {
        //   cout << "dpotrf failed" << endl;
        // }
        // for (i = 0; i < nTime; ++i) {
        //   detB += 2 * log(A[i*nTime + i]);
        // }
        // detB *= nLocations;
        // dpotri_(&lower, &nTime, A, &nTime, &info);
        // if (info) {
        //   cout << "dpotri failed" << endl;
        // }
        // fillUTri(A, nTime);
        // kron(Il, nLocations, nLocations, A, nTime, nTime, B, nTimeLoc, nTimeLoc);
        // fillUTri(lambda, nDay);
        // kron(lambda, nDay, nDay, B, nTimeLoc, nTimeLoc, C, nTimeLocDay, nTimeLocDay);
        // detC = nTimeLoc * detL + nDay * detB;
        // dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
        // logPostThetaCand += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
        
        /************
        * Likelihood
        ***********/
        dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
               eta, &inc);
        // daxpy_(&n, &one, w, &inc, eta, &inc);
        dgemv_(&transN, &n, &nPhi, &one, varDesign, &n, phi, &inc, &zero, 
               sigma, &inc);
        logPostThetaCand += betaLogpostJacobianVar(n, y, eta, sigma);
        
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
      
      // logPostWCurr = 0.0;
      // for (i = 0; i < nPhi; ++i) {
      //   phi[i] = exp(theta[phiIndx + i]);
      // }
      // covTransInvExpand(&theta[lambdaIndx], lambda, nDay);
      // rho = logitInv(theta[rhoIndx], rhoA, rhoB);
      // 
      // detL = 0.0;
      // for (i = 0; i < nDay; ++i) {
      //   detL += 2 * log(lambda[i*nDay + i]);
      // }
      // dpotri_(&lower, &nDay, lambda, &nDay, &info);
      // if(info != 0) {
      //   cout << "c++ error: dpotrf failed" << endl;
      // }
      // AR1(nTime, rho, one, A);
      // detB = 0.0;
      // dpotrf_(&lower, &nTime, A, &nTime, &info);
      // if (info) {
      //   cout << "dpotrf failed" << endl;
      // }
      // for (i = 0; i < nTime; ++i) {
      //   detB += 2 * log(A[i*nTime + i]);
      // }
      // detB *= nLocations;
      // dpotri_(&lower, &nTime, A, &nTime, &info);
      // if (info) {
      //   cout << "dpotri failed" << endl;
      // }
      // fillUTri(A, nTime);
      // kron(Il, nLocations, nLocations, A, nTime, nTime, B, nTimeLoc, nTimeLoc);
      // fillUTri(lambda, nDay);
      // kron(lambda, nDay, nDay, B, nTimeLoc, nTimeLoc, C, nTimeLocDay, nTimeLocDay);
      // detC = nTimeLoc * detL + nDay * detB;
      // dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
      // logPostWCurr += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
      // 
      // /************
      // * Likelihood
      // ***********/
      // dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
      //        eta, &inc);
      // dcopy_(&n, eta, &inc, tmp, &inc);
      // daxpy_(&n, &one, w, &inc, tmp, &inc);
      // dgemv_(&transN, &n, &nPhi, &one, varDesign, &n, phi, &inc, &zero,
      //        sigma, &inc);
      // logPostWCurr += betaLogpostJacobianVar(n, y, tmp, sigma);
      // 
      // /************************************
      // * Proposed random effects
      // ***********************************/
      // for (c = 0; c < n; ++c) {
      //   wCurr = w[c];
      //   w[c] = rnorm(wCurr, exp(wTuning[c]));
      // 
      //   logPostWCand = 0.0;
      // 
      //   dsymv_(&lower, &n, &one, C, &n, w, &inc, &zero, tmp, &inc);
      //   logPostWCand += -0.5 * detC - 0.5 * ddot_(&n, w, &inc, tmp, &inc);
      // 
      //   /************
      //   * Likelihood
      //   ***********/
      //   dcopy_(&n, eta, &inc, tmp, &inc);
      //   daxpy_(&n, &one, w, &inc, tmp, &inc);
      //   logPostWCand += betaLogpostJacobianVar(n, y, tmp, sigma);
      // 
      //   if (runif(0.0, 1.0) <= exp(logPostWCand - logPostWCurr)) {
      //     logPostWCurr = logPostWCand;
      //     wAccept[c]++;
      //   } else {
      //     w[c] = wCurr;
      //   }
      // }
      
      /************************************
      * Save current values
      ***********************************/
      dcopy_(&nTheta, theta, &inc, &thetaSamples[nTheta * d], &inc);
      dcopy_(&n, w ,&inc, &wSamples[n * d], &inc);
      
      dgemv_(&transN, &n, &nSpline, &one, Z, &n, &theta[betaIndx], &inc, &zero,
             eta, &inc);
      // daxpy_(&n, &one, w, &inc, eta, &inc);
      
      for (i = 0; i < nPhi; ++i) {
        phi[i] = exp(theta[phiIndx + i]);
      }
      dgemv_(&transN, &n, &nPhi, &one, varDesign, &n, phi, &inc, &zero, 
             sigma, &inc);
      for (i = 0; i < n; ++i) {
        mu = logitInv(eta[i], 0.0, 1.0);
        precision = sigma[i];
        fittedSamples[d * n + i] = rbeta(mu * precision, (1.0 - mu) * precision);
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
    for (j = 0; j < nPhi; ++j) {
      thetaSamples[i*nTheta + phiIndx + j] = exp(thetaSamples[i*nTheta + phiIndx + j]);
    }
    // covTransInv(&thetaSamples[nTheta*i + lambdaIndx], &thetaSamples[nTheta*i + lambdaIndx],
    //             nDay);
    // thetaSamples[i*nTheta + rhoIndx] = logitInv(thetaSamples[i*nTheta + rhoIndx],
    //                                             rhoA, rhoB);
  }
  
  writeRMatrix(outFile + "-fitted-" + nChain, fittedSamples, n, nSamples);
  writeRMatrix(outFile + "-theta-" + nChain, thetaSamples, nTheta, nSamples);
  writeRMatrix(outFile + "-w-" + nChain, wSamples, n, nSamples);
  
  return 0;
}
