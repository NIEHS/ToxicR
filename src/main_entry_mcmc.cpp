/*
 * Copyright 2020  US. Department of Health and Human Services (HHS),
 * National Institute of Environmental Health Sciences (NIEHS)
 * Email: Matt Wheeler  <matt.wheeler@nih.gov>
 *
 *
 *Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software *and associated documentation files (the "Software"), to deal
 in the Software without restriction, *including without limitation the rights
 to use, copy, modify, merge, publish, distribute, *sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software *is
 furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all
 copies *or substantial portions of the Software.

 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, *INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A *PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF *CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE *OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 *
 *
 */

#include <iostream>
#include <string>
#include <vector>

// const Map<MatrixXd> A(as<Map<MatrixXd>>(AA));

#define STRICT_R_HEADERS

#include <DichGammaBMD_NC.h>
#include <DichHillBMD_NC.h>
#include <DichLogLogisticBMD_NC.h>
#include <DichLogProbitBMD_NC.h>
#include <DichLogisticBMD_NC.h>
#include <DichMultistageBMD_NC.h>
#include <DichProbitBMD_NC.h>
#include <DichQlinearBMD_NC.h>
#include <DichWeibullBMD_NC.h>

#include "mcmc_analysis.h"

#include <bmd_calculate.h>

#include "normal_EXP_NC.h"
#include "normal_HILL_NC.h"
#include "normal_POLYNOMIAL_NC.h"
#include "normal_POWER_NC.h"

#include "lognormal_EXP_NC.h"
#include "lognormal_HILL_NC.h"
#include "lognormal_POLYNOMIAL_NC.h"
#include "lognormal_POWER_NC.h"
#include "seeder.h"

#include "continuous_clean_aux.h"

#ifdef R_COMPILATION
// necessary things to run in R
#ifdef ToxicR_DEBUG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <RcppEigen.h>
#pragma GCC diagnostic pop
#else
#include <RcppEigen.h>
#endif
#include <RcppGSL.h>
#else
#include <Eigen/Dense>
#endif

#include "bmdStruct.h"
#include "bmds_entry.h"
#include "continuous_clean_aux.h"
#include "continuous_entry_code.h"
#include "dichotomous_entry_code.h"

#include "list_r_conversion.h"

using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

Eigen::MatrixXd fix_sample(Eigen::MatrixXd A, dich_model mtype, double max) {

  // Note: Samples are by column.
  switch (mtype) {
  case dich_model::d_hill:
    A.row(2).array() += A.row(3).array() * log(1. / max);
    break;
  case dich_model::d_gamma:
    A.row(2) *= 1. / max;
    break;
  case dich_model::d_logistic:
    A.row(1) *= 1. / max;
    break;
  case dich_model::d_loglogistic:
    A.row(1).array() += A.row(2).array() * log(1. / max);
    break;
  case dich_model::d_logprobit:
    A.row(1).array() += A.row(2).array() * log(1. / max);
    break;
  case dich_model::d_multistage:
    for (int j = 1; j < A.rows(); j++) {
      A.row(j) *= pow(1. / max, j);
    }
    break;
  case dich_model::d_probit:
    A.row(1) *= 1. / max;
    break;
  case dich_model::d_qlinear:
    A.row(1) *= 1. / max;
    break;
  default:
    for (int i = 0; i < A.cols(); i++) {
      A(2, i) *= pow(1 / max, A(1, i));
    }
  }

  return A;
}

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_dichotomous_single_mcmc
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis. Does MCMC sample
// output: BMD analysis with the model specified by NumericVector model
// [[Rcpp::export(".run_dichotomous_single_mcmc")]]
List run_dichotomous_single_mcmc(NumericVector model, Eigen::MatrixXd Y,
                                 Eigen::MatrixXd D, Eigen::MatrixXd pr,
                                 NumericVector options, int seed) {
  Seeder *seeder = Seeder::getInstance();
  seeder->setSeed(seed);
  dichotomous_analysis mcmcAnal;
  mcmcAnal.BMD_type = eExtraRisk; // (options[0]==1)?eExtraRisk:eAddedRisk;
  mcmcAnal.BMR = options[0];
  mcmcAnal.alpha = options[1];
  mcmcAnal.samples = options[2];
  mcmcAnal.burnin = options[3];
  mcmcAnal.parms = pr.rows();
  mcmcAnal.model = (dich_model)model[0];
  mcmcAnal.Y = new double[Y.rows()];
  mcmcAnal.n_group = new double[Y.rows()];
  mcmcAnal.doses = new double[D.rows()];
  mcmcAnal.prior = new double[pr.cols() * pr.rows()];
  mcmcAnal.prior_cols = pr.cols();
  mcmcAnal.n = Y.rows();
  mcmcAnal.degree = 0;

  if (mcmcAnal.model == dich_model::d_multistage) {
    mcmcAnal.degree = mcmcAnal.parms - 1;
  }

  bmd_analysis_MCMC output;
  output.samples = mcmcAnal.samples; // initialize
  output.model = (dich_model)0;
  output.BMDS = new double[mcmcAnal.samples];
  output.parms = new double[mcmcAnal.samples * pr.rows()];

  for (int i = 0; i < Y.rows(); i++) {
    mcmcAnal.Y[i] = Y(i, 0);
    mcmcAnal.n_group[i] = Y(i, 1);
  }

  for (int i = 0; i < D.rows(); i++) {
    mcmcAnal.doses[i] = D(i, 0);
  }

  // copy in column major order
  for (int i = 0; i < pr.rows(); i++) {
    for (int j = 0; j < pr.cols(); j++) {
      mcmcAnal.prior[i + j * pr.rows()] = pr(i, j);
    }
  }

  dichotomous_model_result res;
  res.parms = new double[pr.rows()];
  res.cov = new double[pr.rows() * pr.rows()];
  res.dist_numE = 200;
  res.bmd_dist = new double[res.dist_numE * 2];

  estimate_sm_mcmc(&mcmcAnal, &res, &output);

  List rV = convert_dichotomous_fit_to_list(&res);
  List t2 = convert_MCMC_fit_to_list(&output);

  List data_out =
      List::create(Named("mcmc_result") = t2, Named("fitted_model") = rV);

  delete[] output.BMDS;
  delete[] output.parms;
  delete[] mcmcAnal.Y;
  delete[] mcmcAnal.n_group;
  delete[] mcmcAnal.doses;
  delete[] mcmcAnal.prior;
  delete[] res.parms;
  delete[] res.cov;
  delete[] res.bmd_dist;
  return data_out;
}

//////////////////////////////////////////////////////////////////////////
// function: run_dichotomous_single_mcmc
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis. Does MCMC sample
// output: BMD analysis with the model specified by NumericVector model
// [[Rcpp::export(".run_continuous_single_mcmc")]]
List run_continuous_single_mcmc(NumericVector model, Eigen::MatrixXd Y,
                                Eigen::MatrixXd D, Eigen::MatrixXd priors,
                                NumericVector options, bool is_logNormal,
                                bool suff_stat, int seed) {
  Seeder *seeder = Seeder::getInstance();
  seeder->setSeed(seed);
  unsigned int samples = (unsigned int)options[7];
  unsigned int burnin = (unsigned int)options[6];
  double tail_p = (double)options[2];
  bool bConstVar = (bool)options[5]; // check if it is constant variance
  bool is_increasing = (bool)options[4];
  double alpha = (double)options[3];
  // double bk_prob = (double)options[2];
  double bmrf = (double)options[1];
  int riskType = (int)options[0];
  int transform = options[8];

  continuous_analysis *mcmcAnal = new continuous_analysis;
  distribution dtype;
  if (is_logNormal) {
    dtype = distribution::log_normal;
  } else {
    if (bConstVar) {
      dtype = distribution::normal;
    } else {
      dtype = distribution::normal_ncv;
    }
  }

  mcmcAnal->model = (cont_model)model[0];
  mcmcAnal->Y = new double[Y.rows()];
  mcmcAnal->n = Y.rows();
  mcmcAnal->n_group = new double[Y.rows()];
  mcmcAnal->sd = new double[Y.rows()];
  mcmcAnal->doses = new double[Y.rows()];
  mcmcAnal->prior = new double[priors.rows() * priors.cols()];
  mcmcAnal->isIncreasing = is_increasing;
  mcmcAnal->disttype = dtype;
  mcmcAnal->prior_cols = priors.cols();
  mcmcAnal->parms = priors.rows();
  mcmcAnal->alpha = alpha;
  mcmcAnal->BMD_type = riskType;
  mcmcAnal->BMR = bmrf;
  mcmcAnal->samples = samples;
  mcmcAnal->burnin = burnin;
  mcmcAnal->tail_prob = tail_p;
  mcmcAnal->transform_dose = transform;
  mcmcAnal->suff_stat = suff_stat;
  mcmcAnal->degree = 0;

  //
  // Check on the polynomial stuff
  //
  if (mcmcAnal->model == cont_model::polynomial) {
    // figure out the degree
    if (mcmcAnal->disttype == distribution::normal) {
      mcmcAnal->degree = mcmcAnal->parms - 2;
    } else if (mcmcAnal->disttype == distribution::normal_ncv) {
      mcmcAnal->degree = mcmcAnal->parms - 3;
    } else {
      // throw an error! can'd do log-normal polynomial
      stop("Polynomial-Log-normal models are not allowed.\n Please choose "
           "normal or normal non-constant variance.");
    }
  }

  bmd_analysis_MCMC *output = new bmd_analysis_MCMC;
  output->parms = new double[samples * mcmcAnal->parms];
  output->BMDS = new double[samples];

  ///////////////////////////////////////////////////////////////////////

  for (int i = 0; i < Y.rows(); i++) {
    mcmcAnal->Y[i] = Y(i, 0);
    mcmcAnal->doses[i] = D(i, 0);
    if (suff_stat) {
      mcmcAnal->n_group[i] = Y(i, 1);
      mcmcAnal->sd[i] = Y(i, 2);
    }
  }

  cp_prior(priors, mcmcAnal->prior);
  ////////////////////////////////////

  continuous_model_result *res =
      new_continuous_model_result(mcmcAnal->model, mcmcAnal->parms, 200);

  estimate_sm_mcmc(mcmcAnal, res, output);

  double v_c, v_nc, v_pow;
  estimate_normal_variance(mcmcAnal, &v_c, &v_nc, &v_pow);

  NumericVector v_inform(3);
  v_inform[0] = v_c;
  v_inform[1] = v_nc;
  v_inform[2] = v_pow;
  List rV = convert_continuous_fit_to_list(res);
  List t2 = convert_MCMC_fit_to_list(output);
  List data_out =
      List::create(Named("mcmc_result") = t2, Named("fitted_model") = rV,
                   Named("varOpt") = v_inform);

  del_mcmc_analysis(output);
  del_continuous_model_result(res);
  // del_continuous_analysis(*mcmcAnal);
  delete[] mcmcAnal->Y;
  mcmcAnal->Y = NULL; //
  delete[] mcmcAnal->n_group;
  mcmcAnal->n_group = NULL;
  delete[] mcmcAnal->sd;
  mcmcAnal->sd = NULL;
  delete[] mcmcAnal->doses;
  mcmcAnal->doses = NULL;
  delete[] mcmcAnal->prior;
  mcmcAnal->prior = NULL;
  delete mcmcAnal;

  return wrap(data_out);
}
