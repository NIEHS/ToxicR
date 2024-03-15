/*
 * Copyright 2020  US HHS, NIEHS 
 * Author Matt Wheeler 
 * e-mail: <matt.wheeler@nih.gov> 
 *
 *Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 *and associated documentation files (the "Software"), to deal in the Software without restriction,
 *including without limitation the rights to use, copy, modify, merge, publish, distribute,
 *sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 *is furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all copies
 *or substantial portions of the Software.

 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 *CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 *OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 */

#define STRICT_R_HEADERS
#pragma once
#ifndef mcmc_bmd_calculateH
#define mcmc_bmd_calculateH
#include <iostream>
//#include "stdafx.h"
#include <chrono>
#include <cmath>
#ifndef WIN32
#include <cfloat>
#endif

#ifdef R_COMPILATION
// necessary things to run in R
#include <RcppEigen.h>
#include <RcppGSL.h>
#else
#include <Eigen/Dense>
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>

#include "seeder.h"
#include "dBMDstatmod.h"
#include "cBMDstatmod.h"
#include "bmd_calculate.h"
#include "mcmc_struct.h"
#include <sys/time.h>
#include <bmdStruct.h>

using namespace std;

#define MCMC_PROPOSAL(MU, X, COV) \
  -0.5 * (MU.transpose() - X.transpose()) * COV.inverse() * (MU - X) // note the integrating constant cancels out

/**********************************************************************
 *  function: MCMC_bmd_analysis_DNC(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
 *							 double BMR, bool isExtra, double alpha, double step_size)
 *  Purpose: provide a MCMC analysis of the specified BMD analysis.
 * *******************************************************************/
template <class LL, class PR>
mcmcSamples MCMC_bmd_analysis_DNC(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
                                  std::vector<bool> fixedB, std::vector<double> fixedV, int degree,
                                  double BMR, bool isExtra, double alpha, int samples)
{

  LL dichotimousM(Y, D, degree);
  PR model_prior(prior);

  dBMDModel<LL, PR> model(dichotimousM, model_prior, fixedB, fixedV);
  optimizationResult oR = findMAP<LL, PR>(&model);

  mcmcSamples rVal;
  rVal.map = oR.functionV;
  rVal.map_estimate = oR.max_parms;
  ;
  rVal.map_cov = model.varMatrix(oR.max_parms);
  rVal.isExtra = false;

  int n = oR.max_parms.rows();
  // generate random univariate normals for the MCMC sampler
  //  there are n x samples generated
  //  ziggurat is used as it is the fastest sampler algorithm gsl has
  Eigen::MatrixXd rNormal(n, samples);
  rNormal.setZero(); 
  Seeder* seeder = Seeder::getInstance();
  for (int i = 0; i < samples; i++)
  {
    for (int j = 0; j < n; j++)
    {
      rNormal(j, i) = seeder->get_gaussian_ziggurat();
    }
  }

  // now sample, samples, number of proposals for the
  // metropolis sampler.
  /*
  FullPivLU<Matrix3f> lu_decomp(your_matrix);
  auto rank = lu_decomp.rank();
  */
  double eps = 1.25;
  Eigen::MatrixXd mu = oR.max_parms;
  Eigen::MatrixXd cov = pow(eps, 2) * model.varMatrix(oR.max_parms);

  Eigen::MatrixXd chol = cov.llt().matrixL();
  Eigen::MatrixXd nSamples = chol * rNormal; // variance of each row
                                             // is is now L'L = cov
  Eigen::MatrixXd penLike(1, samples);
  penLike.setZero(); 
  Eigen::MatrixXd BMD(1, samples);
  BMD.setZero(); 
  /////////////////////////////////////////////////////////////////
  nSamples.col(0) = mu; // starting value of the MCMC

  penLike(0, 0) = -model.negPenLike(nSamples.col(0));
  for (int i = 1; i < samples; i++)
  {
    Eigen::MatrixXd value = nSamples.col(i) + nSamples.col(i - 1); // mu;
    // Metropolis
    //  make sure the proposal isn't imposible

    double t = model.prior_model.neg_log_prior(value);

    if (!isnan(t) &&
        !isinf(t))
    {
      Eigen::MatrixXd a = MCMC_PROPOSAL(value, nSamples.col(i - 1), cov);
      Eigen::MatrixXd b = MCMC_PROPOSAL(nSamples.col(i - 1), value, cov);
      double numer = -model.negPenLike(value) + a(0, 0);
      double denom = penLike(0, i - 1) + b(0, 0);

      double test = exp(numer - denom);
      double rr = seeder->get_uniform();

      if (isnan(test))
      {
        test = 0.0; // no probability of making this jump
      }

      if (rr < test)
      {
        nSamples.col(i) = value;
        penLike(0, i) = -model.negPenLike(value);
      }
      else
      {
        nSamples.col(i) = nSamples.col(i - 1);
        penLike(0, i) = penLike(0, i - 1);
      }
      // compute the calculated BMD for each sample
    }
    else
    {
      // previous proposal was outside of the bounds
      nSamples.col(i) = nSamples.col(i - 1);
      penLike(0, i) = penLike(0, i - 1);
    }
    BMD(0, i) = isExtra ? model.log_likelihood.compute_BMD_EXTRA_NC(nSamples.col(i), BMR)
                        : model.log_likelihood.compute_BMD_ADDED_NC(nSamples.col(i), BMR);
  }

  /////////////////////////////////////////////////////////////////
  rVal.BMD = BMD;          // vector of burn in BMD samples
  rVal.samples = nSamples; // vector of sample parameters
  rVal.isExtra = isExtra;
  rVal.log_posterior = penLike;
  rVal.BMR = BMR;
  //////////////////////////////////////////////////////////////////

  return rVal;
}

/**********************************************************************
 *  function: MCMC_bmd_analysis_DNC(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
 *							 double BMR, bool isExtra, double alpha, double step_size)
 *  Purpose: provide a MCMC analysis of the specified BMD analysis.
 ********************************************************************/
template <class LL, class PR>
mcmcSamples mcmc_continuous(cBMDModel<LL, PR> *model, int samples,
                            contbmd BMDt, double BMR, double pointP,
                            Eigen::MatrixXd init_vals = Eigen::MatrixXd::Zero(1, 1))
{

  optimizationResult oR;
  if (init_vals.rows() == 1 &&
      init_vals.cols() == 1 &&
      init_vals(0, 0) == 0.0)
  {
    oR = findMAP<LL, PR>(model);
  }
  else if (model->modelling_type() == cont_model::gamma_aerts)
  {
	oR = findMAP<LL, PR>(model);
  }
  else
  {
    oR = findMAP<LL, PR>(model, init_vals);
  }

  mcmcSamples rVal;
  rVal.isExtra = false;
  rVal.map = oR.functionV;
  rVal.map_estimate = oR.max_parms;
  ;
//  cout << "MCMC findMAP init = " << init_vals << endl;

  rVal.map_cov = model->varMatrix(oR.max_parms);

//  cout << "MCMC cov";

  struct timeval tv; // Seed generation based on time
  gettimeofday(&tv, 0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;

  int n = oR.max_parms.rows();
  // generate random univariate normals for the MCMC sampler
  //  there are n x samples generated
  //  ziggurat is used as it is the fastest sampler algorithm gsl has
  Eigen::MatrixXd rNormal(n, samples);
  rNormal.setZero(); 

  Seeder* seeder = Seeder::getInstance();
  for (int i = 0; i < samples; i++)
  {
    for (int j = 0; j < n; j++)
    {
      rNormal(j, i) = seeder->get_gaussian_ziggurat();
    }
  }

  // now sample From a metropolis-Hastings sampler.
  double eps = 1.25;
  Eigen::MatrixXd mu = oR.max_parms;
//  cout << "MCMC mean " << mu;
  Eigen::MatrixXd cov = pow(eps, 2) * model->varMatrix(oR.max_parms);
  // if there is a 0 column, add value so cholesky is computed
  for(unsigned i = 0; i < n; i++){
    if(cov(i,i) == 0.0){
      cov(i,i) += 1e6;
    }
  }
  Eigen::MatrixXd chol = cov.llt().matrixL();
  //undo the cholesky addition
  for(unsigned i = 0; i < n; i++){
    if(chol(i,i) == 1e3){
      chol(i,i) = 0.0;
    }
  }
  Eigen::MatrixXd nSamples = chol * rNormal; // variance of each row
  // is is now LL' = cov
  Eigen::MatrixXd penLike(1, samples);
  penLike.setZero(); 
  Eigen::MatrixXd BMD(1, samples);
  BMD.setZero(); 
  /////////////////////////////////////////////////////////////////
  nSamples.col(0) = mu; // starting value of the MCMC

  penLike(0, 0) = -model->negPenLike(nSamples.col(0));

  for (int i = 1; i < samples; i++)
  { // nSamples.col(i-1);
    // Metropolis-Hasings proposal
    //  make sure the proposal isn't imposible
    Eigen::MatrixXd value = nSamples.col(i) + nSamples.col(i - 1);
    double t = model->prior_model.neg_log_prior(value);

    if (!isnan(t) &&
        !isinf(t))
    {
      Eigen::MatrixXd a = MCMC_PROPOSAL(value, nSamples.col(i - 1), cov);
      Eigen::MatrixXd b = MCMC_PROPOSAL(nSamples.col(i - 1), value, cov);
      double numer = -model->negPenLike(value) + a(0, 0);
      double denom = penLike(0, i - 1) + b(0, 0);
      double test = exp(numer - denom);
      double rr = seeder->get_uniform();

      if (isnan(test))
      {
        test = 0.0; // no probability of making this jump
                    // cout << i << endl;
      }

      if (rr < test)
      {
        nSamples.col(i) = value;
        penLike(0, i) = -model->negPenLike(value);
      }
      else
      {
        nSamples.col(i) = nSamples.col(i - 1);
        penLike(0, i) = penLike(0, i - 1);
      }
    }
    else
    {
      // previous proposal was outside of the bounds
      nSamples.col(i) = nSamples.col(i - 1);
      penLike(0, i) = penLike(0, i - 1);
    }

    BMD(0, i) = model->returnBMD(nSamples.col(i), BMDt,
                                 BMR, pointP);
  }

  /////////////////////////////////////////////////////////////////
  rVal.BMD = BMD;          // vector of BMD samples
  rVal.samples = nSamples; // vector of sample parameters
  rVal.log_posterior = penLike;
  rVal.BMR = BMR;
  //////////////////////////////////////////////////////////////////
  return rVal;
}

template <class LL, class PR>
mcmcSamples MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
                                                   std::vector<bool> fixedB, std::vector<double> fixedV,
                                                   bool isIncreasing, double point_p, bool suff_stat,
                                                   double BMR, contbmd BMDType, double alpha,
                                                   int samples, int adverse_d,
                                                   Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1, 1))
{

  LL likelihood(Y, D, suff_stat, adverse_d);
  // cout << prior << endl;
  PR model_prior(prior);

  cBMDModel<LL, PR> model(likelihood, model_prior, fixedB, fixedV, isIncreasing);
  return mcmc_continuous<LL, PR>(&model, samples, BMDType, BMR, point_p, init);
}

template <class LL, class PR>
mcmcSamples MCMC_bmd_analysis_CONTINUOUS_NORMAL(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
                                                std::vector<bool> fixedB, std::vector<double> fixedV,
                                                bool isIncreasing, double point_p, bool suff_stat,
                                                double BMR, contbmd BMDType, bool const_var,
                                                double alpha, int samples, int adverse_d,
                                                Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1, 1))
{

  LL likelihood(Y, D, suff_stat, const_var, adverse_d);

  PR model_prior(prior);

  cBMDModel<LL, PR> model(likelihood, model_prior, fixedB, fixedV, isIncreasing);

  return mcmc_continuous<LL, PR>(&model, samples, BMDType, BMR, point_p, init);
}

#endif
