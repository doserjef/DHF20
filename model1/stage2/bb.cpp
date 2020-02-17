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
  cout << nChain << endl;
  int nThreads; par.getVal("n.threads", nThreads);
  string outFile; par.getVal("outfile", outFile);
  
  int nLocations; par.getVal("n.locations", nLocations);
  int nTime; par.getVal("n.time", nTime);
  int nSpline; par.getVal("n.spline", nSpline);
  int nPhi; par.getVal("n.phi", nPhi);
  
  double *x; int n, a; par.getFile("x.file", x, n, a);
  double *y; int colY; par.getFile("y.file", y, n, colY);
  
  vector<double> betaTuning; par.getVal("beta.tuning", betaTuning);
  vector<double> betaStarting; par.getVal("beta.starting", betaStarting);
  vector<double> betaMu; par.getVal("beta.mu", betaMu);
  vector<double> betaSd; par.getVal("beta.sd", betaSd);
  
  double phiTuning; par.getVal("phi.tuning", phiTuning);
  double phiStarting; par.getVal("phi.starting", phiStarting);
  double phiA; par.getVal("phi.a", phiA);
  double phiB; par.getVal("phi.b", phiB);
  
  double sigmaSqTuning; par.getVal("sigma.sq.tuning", sigmaSqTuning);
  double sigmaSqStarting; par.getVal("sigma.sq.starting", sigmaSqStarting);
  double sigmaSqA; par.getVal("sigma.sq.a", sigmaSqA);
  double sigmaSqB; par.getVal("sigma.sq.b", sigmaSqB);
  
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
  char lower, transN;
  int info, inc;
  int nTimeSq, nBreak;
  int nTheta, betaIndx, sigmaSqIndx, rhoIndx, phiIndx; 
  double alphaMax, alphaMin, mu, thetaCurr, logPostWCurr;
  double logPostThetaCurr, logPostThetaCand, sigmaSq, rho, logDet, zero, one;
  double logPostWCand, wCurr, acceptRate, phi;
  double *Z, *A, *w, *tmp, *wTuning, *eta, *wAccept;
  double *theta, *thetaTuning, *accept;
  double *thetaSamples, *wSamples, *fittedSamples, *zSamples;
  
  /****************************************************************************
   * Log tuning parameters for use in adaptive MCMC sampler
   ***************************************************************************/
  nTheta = nSpline + nPhi + 1 + 1;
  betaIndx = 0;
  thetaTuning = new double[nTheta];
  phiIndx = betaIndx + nSpline;
  sigmaSqIndx = phiIndx + nPhi;
  rhoIndx = sigmaSqIndx + 1;
  for (i =0; i < nSpline; ++i) {
    thetaTuning[betaIndx + i] = log(betaTuning[i]);
  }
  thetaTuning[phiIndx] = log(phiTuning);
  thetaTuning[sigmaSqIndx] = log(sigmaSqTuning);
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
  theta[sigmaSqIndx] = log(sigmaSqStarting);
  theta[rhoIndx] = logit(rhoStarting, rhoA, rhoB);
  
  nTimeSq = nTime * nTime;
  A = new double[nTimeSq];
  lower = 'L';
  info = 0;
  w = new double[n]; zeros(w, n);
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
  
  /****************************************************************************
   * Spline Coefficients
   ***************************************************************************/
  gsl_bspline_workspace *bw; 
  nBreak = nSpline - 2;
  bw = gsl_bspline_alloc(4, nBreak);
  gsl_vector *B = gsl_vector_alloc(nSpline);
  alphaMin = 1.0;
  alphaMax = 0.0;
  gsl_vector *knots = gsl_vector_alloc(nBreak);
  gsl_vector_set(knots, 0, 0.0);
  gsl_vector_set(knots, 1, 0.7);
  gsl_vector_set(knots, 2, 0.8);
  gsl_vector_set(knots, 3, 0.97);
  gsl_vector_set(knots, 4, 0.99);
  gsl_vector_set(knots, 5, 1.0);
  // // Manually set breakpoints
  gsl_bspline_knots(knots, bw);
  Z = new double[n * nSpline];
  
  // for (i = 0; i < n; ++i) {
  //   gsl_bspline_eval(x[i], B, bw);
  //   for (j = 0; j < nSpline; ++j) {
  //     Z[j * n + i] = gsl_vector_get(B, j);
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
        gsl_bspline_eval(x[d*n + i], B, bw);
        for (j = 0; j < nSpline; ++j) {
          Z[j * n + i] = gsl_vector_get(B, j);
        }
      }
      /************************************
       * Current theta
       ***********************************/
      logPostThetaCurr = 0.0;
      phi = exp(theta[phiIndx]);
      sigmaSq = exp(theta[sigmaSqIndx]);
      rho = logitInv(theta[rhoIndx], rhoA, rhoB);

      /************
       * Priors
       ***********/
      for (i = 0; i < nSpline; ++i) {
        logPostThetaCurr += dnorm(theta[betaIndx + i], betaMu[i], betaSd[i], 1);
      }
      logPostThetaCurr += iGammaLogpostJacobian(phi, phiA, phiB);
      logPostThetaCurr += iGammaLogpostJacobian(sigmaSq, sigmaSqA, sigmaSqB);
      logPostThetaCurr += uniformLogpostJacobian(rho, rhoA, rhoB);

      AR1(nTime, rho, sigmaSq, A);
      logDet = 0.0;
      dpotrf_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotrf failed" << endl;
      }
      for (i = 0; i < nTime; ++i) {
        logDet += 2 * log(A[i*nTime + i]);
      }
      logDet *= nLocations;
      dpotri_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotri failed" << endl;
      }

      for (i = 0; i < nLocations; ++i) {
        dsymv_(&lower, &nTime, &one, A, &nTime, &w[nTime * i], &inc, &zero,
               &tmp[nTime * i], &inc);
      }
      logPostThetaCurr += -0.5 * logDet - 0.5 * ddot_(&n, w, &inc, tmp, &inc);

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
        sigmaSq = exp(theta[sigmaSqIndx]);
        rho = logitInv(theta[rhoIndx], rhoA, rhoB);

        /************
         * Priors
         ***********/
        for (i = 0; i < nSpline; ++i) {
          logPostThetaCand += dnorm(theta[betaIndx + i], betaMu[i], betaSd[i], 1);
        }
        logPostThetaCand += iGammaLogpostJacobian(phi, phiA, phiB);
        logPostThetaCand += iGammaLogpostJacobian(sigmaSq, sigmaSqA, sigmaSqB);
        logPostThetaCand += uniformLogpostJacobian(rho, rhoA, rhoB);

        AR1(nTime, rho, sigmaSq, A);
        logDet = 0.0;
        dpotrf_(&lower, &nTime, A, &nTime, &info);
        if (info) {
          cout << "dpotrf failed" << endl;
        }
        for (i = 0; i < nTime; ++i) {
          logDet += 2 * log(A[i*nTime + i]);
        }
        logDet *= nLocations;
        dpotri_(&lower, &nTime, A, &nTime, &info);
        if (info) {
          cout << "dpotri failed" << endl;
        }

        for (i = 0; i < nLocations; ++i) {
          dsymv_(&lower, &nTime, &one, A, &nTime, &w[nTime * i], &inc, &zero,
                 &tmp[nTime * i], &inc);
        }
        logPostThetaCand += -0.5 * logDet - 0.5 * ddot_(&n, w, &inc, tmp, &inc);

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
      sigmaSq = exp(theta[sigmaSqIndx]);
      rho = logitInv(theta[rhoIndx], rhoA, rhoB);
      AR1(nTime, rho, sigmaSq, A);
      logDet = 0.0;
      dpotrf_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotrf failed" << endl;
      }
      for (i = 0; i < nTime; ++i) {
        logDet += 2 * log(A[i*nTime + i]);
      }
      logDet *= nLocations;
      dpotri_(&lower, &nTime, A, &nTime, &info);
      if (info) {
        cout << "dpotri failed" << endl;
      }

      for (i = 0; i < nLocations; ++i) {
        dsymv_(&lower, &nTime, &one, A, &nTime, &w[nTime * i], &inc, &zero,
               &tmp[nTime * i], &inc);
      }
      logPostWCurr += -0.5 * logDet - 0.5 * ddot_(&n, w, &inc, tmp, &inc);

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

        for (i = 0; i < nLocations; ++i) {
          dsymv_(&lower, &nTime, &one, A, &nTime, &w[nTime * i], &inc, &zero,
                 &tmp[nTime * i], &inc);
        }
        logPostWCand += -0.5 * logDet - 0.5 * ddot_(&n, w, &inc, tmp, &inc);

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
    thetaSamples[i*nTheta + sigmaSqIndx] = exp(thetaSamples[i*nTheta + sigmaSqIndx]);
    thetaSamples[i*nTheta + rhoIndx] = logitInv(thetaSamples[i*nTheta + rhoIndx],
                                                rhoA, rhoB);
  }

  writeRMatrix(outFile + "-fitted-" + nChain, fittedSamples, n, nSamples);
  writeRMatrix(outFile + "-theta-" + nChain, thetaSamples, nTheta, nSamples);
  writeRMatrix(outFile + "-w-" + nChain, wSamples, n, nSamples);
  
  return 0;
}
