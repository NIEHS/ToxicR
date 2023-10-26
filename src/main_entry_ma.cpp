/*
 * Copyright 2020  US. Department of Health and Human Services (HHS), 
 * National Institute of Environmental Health Sciences (NIEHS)
 * Email: Matt Wheeler  <matt.wheeler@nih.gov>
 *
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
#include <RcppEigen.h>
#include <RcppGSL.h>

#include <string>
#include <vector>
#include <limits>
#include <cmath>
//#include <math.h>
#include <stdio.h>
#include "bmds_entry.h"
#include "bmdStruct.h"
#include "continuous_model_functions.h"

using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

// const Map<MatrixXd> A(as<Map<MatrixXd>>(AA));

#include <statmod.h>

#include <log_likelihoods.h>
#include <normal_likelihoods.h>
#include <normalModels.h>
#include <binomModels.h>
#include <IDPrior.h>

#include <bmd_calculate.h>

#include "normal_HILL_NC.h"
#include "normal_POWER_NC.h"
#include "normal_POLYNOMIAL_NC.h"
#include "normal_EXP_NC.h"

#include "lognormal_HILL_NC.h"
#include "lognormal_POWER_NC.h"
#include "lognormal_POLYNOMIAL_NC.h"
#include "lognormal_EXP_NC.h"

#include "continuous_clean_aux.h"
#include "continuous_entry_code.h"
#include "dichotomous_entry_code.h"
/*
 *
 */
List convert_dichotomous_fit_to_list(dichotomous_model_result *result)
{
  NumericVector parms(result->nparms);
  NumericMatrix covM(result->nparms, result->nparms);

  for (int i = 0; i < result->nparms; i++)
  {
    parms[i] = result->parms[i];
    for (int j = 0; j < result->nparms; j++)
    {
      covM(i, j) = result->cov[i + j * result->nparms];
    }
  }
  char str[160];

  switch (result->model)
  {

  case dich_model::d_hill:
    snprintf(str,160, "Model:  %s", "Hill");
    break;
  case dich_model::d_gamma:
    snprintf(str,160, "Model:  %s", "Gamma");
    break;
  case dich_model::d_logistic:
    snprintf(str,160, "Model:  %s", "Logistic");
    break;
  case dich_model::d_loglogistic:
    snprintf(str,160, "Model:  %s", "Log-Logistic");
    break;
  case dich_model::d_logprobit:
    snprintf(str,160, "Model:  %s", "Log-Probit");
    break;
  case dich_model::d_multistage:
    snprintf(str,160, "Model:  %s", "Multistage");
    break;
  case dich_model::d_qlinear:
    snprintf(str,160, "Model:  %s", "Quantal-Linear");
    break;
  case dich_model::d_probit:
    snprintf(str,160, "Model:  %s", "Probit");
    break;
  case dich_model::d_weibull:
    snprintf(str,160, "Model: %s", "Weibull");
    break;
  default:
    snprintf(str,160, "Model:  %s", "Danger");
    break;
  }
  double maximum = result->max;
  NumericMatrix bmd_distribution(result->dist_numE, 2);

  for (int i = 0; i < result->dist_numE; i++)
  {
    bmd_distribution(i, 0) = result->bmd_dist[i];
    bmd_distribution(i, 1) = result->bmd_dist[i + result->dist_numE];
  }

  List rV = List::create(Named("full_model") = str,
                         Named("parameters") = parms,
                         Named("covariance") = covM,
                         Named("bmd_dist") = bmd_distribution,
                         Named("bmd") = result->bmd,
                         Named("maximum") = maximum,
                         Named("gof_p_value") = result->gof_p_value,
                         Named("gof_chi_sqr_statistic") = result->gof_chi_sqr_statistic);

  rV.attr("class") = "BMDdich_fit_maximized";
  return rV;
}

List convert_dichotomous_maresults_to_list(dichotomousMA_result *result)
{
  List fittedModels;
  char str[160];

  for (int i = 0; i < result->nmodels; i++)
  {
    snprintf(str,160, "Fitted_Model_%d", i + 1);
    fittedModels.push_back(convert_dichotomous_fit_to_list(result->models[i]),
                           str);
  }

  NumericMatrix ma_bmd_dist(result->dist_numE, 2);
  NumericVector post_probs(result->nmodels);
  for (int i = 0; i < result->dist_numE; i++)
  {
    ma_bmd_dist(i, 0) = result->bmd_dist[i];
    ma_bmd_dist(i, 1) = result->bmd_dist[i + result->dist_numE];
  }
  for (int i = 0; i < result->nmodels; i++)
  {
    post_probs[i] = result->post_probs[i];
  }

  fittedModels.push_back(ma_bmd_dist, "BMD_CDF");
  fittedModels.push_back(post_probs, "posterior_probs");

  return fittedModels;
}

/*
 *
 *
 */
List convert_continuous_fit_to_list(continuous_model_result *result)
{
  NumericVector parms(result->nparms);
  NumericMatrix covM(result->nparms, result->nparms);

  for (int i = 0; i < result->nparms; i++)
  {
    parms[i] = result->parms[i];
    for (int j = 0; j < result->nparms; j++)
    {
      covM(i, j) = result->cov[i + j * result->nparms];
    }
  }
  char dist[160];
  char str[360];

  switch (result->dist)
  {
  case distribution::normal:
    snprintf(dist,160, "Distribution: %s", "Normal");
    break;
  case distribution::normal_ncv:
    snprintf(dist,160, "Distribution: %s", "Normal-NCV");
    break;
  case distribution::log_normal:
    snprintf(dist,160, "Distribution: %s", "Log-Normal");
    break;
  }

  switch (result->model)
  {
  case cont_model::hill:
    snprintf(str,360, "Model: %s %s", "Hill", dist);
    break;
  case cont_model::exp_3:
    snprintf(str,360, "Model: %s %s", "Exponential-3", dist);
    break;
  case cont_model::exp_5:
    snprintf(str,360, "Model: %s %s", "Exponential-5", dist);
    break;
  case cont_model::power:
    snprintf(str,360, "Model: %s %s", "Power", dist);
    break;
  case cont_model::funl:
    snprintf(str,360, "Model: %s %s", "FUNL", dist);
    break;
  case cont_model::polynomial:
    snprintf(str,360, "Model: %s %s", "Polynomial", dist);
    break;
  case cont_model::exp_aerts:
    snprintf(str,360, "Model: %s %s", "Exponential-Aerts", dist);
    break;
  case cont_model::invexp_aerts:
    snprintf(str,360, "Model: %s %s", "Inverse Exponential-Aerts", dist);
    break;
  case cont_model::gamma_aerts:
    snprintf(str,360, "Model: %s %s", "Gamma-Aerts", dist);
    break;
  case cont_model::invgamma_aerts:
    snprintf(str,360, "Model: %s %s", "Inverse Gamma-Aerts", dist);
    break;
  case cont_model::hill_aerts:
    snprintf(str,360, "Model: %s %s", "Hill-Aerts", dist);
    break;
  case cont_model::lomax_aerts:
    snprintf(str,360, "Model: %s %s", "Lomax-Aerts", dist);
    break;
  case cont_model::invlomax_aerts:
    snprintf(str,360, "Model: %s %s", "Inverse Lomax-Aerts", dist);
    break;
  case cont_model::lognormal_aerts:
    snprintf(str,360, "Model: %s %s", "Lognormal-Aerts", dist);
    break;
  case cont_model::logskew_aerts:
    snprintf(str,360, "Model: %s %s", "Log-skew-normal-Aerts", dist);
    break;
  case cont_model::invlogskew_aerts:
    snprintf(str,360, "Model: %s %s", "Inverse Log-skew-normal-Aerts", dist);
    break;
  case cont_model::logistic_aerts:
    snprintf(str,360, "Model: %s %s", "Logistic-Aerts", dist);
    break;
  case cont_model::probit_aerts:
    snprintf(str,360, "Model: %s %s", "Probit-Aerts", dist);
    break;
  case cont_model::LMS:
    snprintf(str,360, "Model: %s %s", "LMS", dist);
    break;
  case cont_model::gamma_efsa:
    snprintf(str,360, "Model: %s %s", "Gamma-EFSA", dist);
    break;
  default:
    snprintf(str,360, "Model: %s %s", "Danger", "Danger");
    break;
  }
  double maximum = result->max;
  NumericMatrix bmd_distribution(result->dist_numE, 2);

  for (int i = 0; i < result->dist_numE; i++)
  {
    bmd_distribution(i, 0) = result->bmd_dist[i];
    bmd_distribution(i, 1) = result->bmd_dist[i + result->dist_numE];
  }

  List rV = List::create(Named("full_model") = str,
                         Named("bmd") = result->bmd,
                         Named("parameters") = parms,
                         Named("covariance") = covM,
                         Named("bmd_dist") = bmd_distribution,
                         Named("maximum") = maximum);
  return rV;
}

/////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////

List convert_continuous_maresults_to_list(continuousMA_result *result)
{

  List fittedModels;
  char str[160];

  for (int i = 0; i < result->nmodels; i++)
  {
    snprintf(str,160, "Fitted_Model_%d", i + 1);
    fittedModels.push_back(convert_continuous_fit_to_list(result->models[i]),
                           str);
  }
  NumericMatrix ma_bmd_dist(result->dist_numE, 2);
  NumericVector post_probs(result->nmodels);
  for (int i = 0; i < result->dist_numE; i++)
  {
    ma_bmd_dist(i, 0) = result->bmd_dist[i];
    ma_bmd_dist(i, 1) = result->bmd_dist[i + result->dist_numE];
  }

  for (int i = 0; i < result->nmodels; i++)
  {
    post_probs[i] = result->post_probs[i];
  }

  fittedModels.push_back(ma_bmd_dist, "ma_bmd");
  fittedModels.push_back(post_probs, "posterior_probs");

  return fittedModels;
}
/////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export(".run_continuous_ma_laplace")]]
List run_continuous_ma_laplace(List model_priors, NumericVector model_type,
                               NumericVector dist_type,
                               Eigen::MatrixXd Y, Eigen::MatrixXd X,
                               NumericVector options)
{

  bool is_increasing = (bool)options[4];
  // double alpha = (double)options[3];
  double tail_p = (double)options[2];
  double bmrf = (double)options[1];
  int riskType = (int)options[0];
  unsigned int samples = (unsigned int)options[5];

  continuousMA_analysis ma_anal;

  ma_anal.nmodels = model_priors.length();
  ma_anal.modelPriors = new double[ma_anal.nmodels];
  ma_anal.priors = new double *[ma_anal.nmodels];
  ma_anal.nparms = new int[ma_anal.nmodels];
  ma_anal.actual_parms = new int[ma_anal.nmodels];
  ma_anal.prior_cols = new int[ma_anal.nmodels];
  ma_anal.models = new int[ma_anal.nmodels];
  ma_anal.disttype = new int[ma_anal.nmodels];

  continuousMA_result *ma_result = new continuousMA_result;
  ma_result->nmodels = ma_anal.nmodels;
  ma_result->dist_numE = 400;
  ma_result->bmd_dist = new double[400 * 2];
  ma_result->post_probs = new double[ma_anal.nmodels];
  ma_result->models = new continuous_model_result *[ma_anal.nmodels];

  for (int i = 0; i < ma_anal.nmodels; i++)
  {
    ma_anal.modelPriors[i] = 1.0 / double(ma_anal.nmodels);
    Eigen::MatrixXd temp = model_priors[i];
    ma_anal.priors[i] = new double[temp.rows() * temp.cols()];

    cp_prior(temp, ma_anal.priors[i]);
    ma_anal.nparms[i] = temp.rows();
    ma_anal.prior_cols[i] = temp.cols();
    ma_anal.models[i] = (int)model_type[i];
    ma_anal.disttype[i] = (int)dist_type[i];
    ma_result->models[i] = new_continuous_model_result(ma_anal.models[i],
                                                       ma_anal.nparms[i],
                                                       400); // have 200 equally spaced values
  }

  /// Set up the other info
  continuous_analysis anal;
  anal.Y = new double[Y.rows()];
  anal.n = Y.rows();
  anal.n_group = new double[Y.rows()];
  anal.sd = new double[Y.rows()];
  anal.doses = new double[Y.rows()];
  anal.prior = NULL;
  anal.isIncreasing = is_increasing;
  anal.alpha = 0.005; // alpha for analyses;
  anal.BMD_type = riskType;
  anal.BMR = bmrf;
  anal.samples = samples;
  anal.tail_prob = tail_p;
  anal.suff_stat = Y.cols() == 3;

  for (int i = 0; i < Y.rows(); i++)
  {
    anal.Y[i] = Y(i, 0);
    anal.doses[i] = X(i, 0);
    if (Y.cols() == 3)
    { // sufficient statistics
      anal.n_group[i] = Y(i, 1);
      anal.sd[i] = Y(i, 2);
    }
  }

  estimate_ma_laplace(&ma_anal, &anal, ma_result);
  List rV = convert_continuous_maresults_to_list(ma_result);

  // free up memory
  for (int i = 0; i < ma_result->nmodels; i++)
  {
    del_continuous_model_result(ma_result->models[i]);
  }
  delete[] ma_result->models;
  delete[] ma_result->post_probs;
  delete[] ma_result->bmd_dist;
  delete ma_result;
  del_continuous_analysis(anal);
  del_continuousMA_analysis(ma_anal);
  return rV;
}

List convert_MCMC_fit_to_list(bmd_analysis_MCMC *a)
{
  List rV;
  NumericMatrix parameters(a->samples, a->nparms);
  NumericMatrix BMDS(a->samples, 1);

  for (unsigned int i = 0; i < a->samples; i++)
  {
    BMDS[i] = a->BMDS[i];
    for (unsigned int j = 0; j < a->nparms; j++)
    {
      parameters(i, j) = a->parms[i + j * a->samples];
    }
  }
  rV = List::create(Named("BMD_samples") = BMDS, Named("PARM_samples") = parameters);
  return rV;
}

List convert_mcmc_results(const ma_MCMCfits *a)
{
  List rV;
  char str[160];

  for (unsigned int i = 0; i < a->nfits; i++)
  {
    snprintf(str,160, "Fitted_Model_%d", i + 1);
    rV.push_back(convert_MCMC_fit_to_list(a->analyses[i]),
                 str);
  }
  return rV;
}
/////////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export(".run_continuous_ma_mcmc")]]
List run_continuous_ma_mcmc(List model_priors, NumericVector model_type,
                            NumericVector dist_type,
                            Eigen::MatrixXd Y, Eigen::MatrixXd X,
                            NumericVector options)
{
  unsigned int burnin = (unsigned int)options[6];
  bool is_increasing = (bool)options[4];
  // double alpha = (double)options[3];
  double tail_p = (double)options[2];
  double bmrf = (double)options[1];
  int riskType = (int)options[0];
  unsigned int samples = (unsigned int)options[5];

  continuousMA_analysis ma_anal;

  ma_anal.nmodels = model_priors.length();
  ma_anal.modelPriors = new double[ma_anal.nmodels];
  ma_anal.priors = new double *[ma_anal.nmodels];
  ma_anal.nparms = new int[ma_anal.nmodels];
  ma_anal.actual_parms = new int[ma_anal.nmodels];
  ma_anal.prior_cols = new int[ma_anal.nmodels];
  ma_anal.models = new int[ma_anal.nmodels];
  ma_anal.disttype = new int[ma_anal.nmodels];

  continuousMA_result *ma_result = new continuousMA_result;
  ma_MCMCfits model_mcmc_info;
  model_mcmc_info.analyses = new bmd_analysis_MCMC *[ma_anal.nmodels];
  model_mcmc_info.nfits = ma_anal.nmodels;
  ma_result->nmodels = ma_anal.nmodels;
  ma_result->dist_numE = 600;
  ma_result->bmd_dist = new double[600 * 2];
  ma_result->post_probs = new double[ma_anal.nmodels];
  ma_result->models = new continuous_model_result *[ma_anal.nmodels];

  for (int i = 0; i < ma_anal.nmodels; i++)
  {
    ma_anal.modelPriors[i] = 1.0 / double(ma_anal.nmodels);
    Eigen::MatrixXd temp = model_priors[i];
    ma_anal.priors[i] = new double[temp.rows() * temp.cols()];
    // cout << temp << endl;
    cp_prior(temp, ma_anal.priors[i]);
    ma_anal.nparms[i] = temp.rows();
    ma_anal.prior_cols[i] = temp.cols();
    ma_anal.models[i] = (int)model_type[i];
    ma_anal.disttype[i] = (int)dist_type[i];
    // cout << ma_anal.models[i] << " " << dist_type[i] << endl;
    ma_result->models[i] = new_continuous_model_result(ma_anal.models[i],
                                                       ma_anal.nparms[i],
                                                       300); // have 300 equally spaced values
    model_mcmc_info.analyses[i] = new_mcmc_analysis(ma_anal.models[i],
                                                    ma_anal.nparms[i],
                                                    samples);
  }

  /// Set up the other info
  continuous_analysis anal;
  anal.Y = new double[Y.rows()];
  anal.n = Y.rows();
  anal.n_group = new double[Y.rows()];
  anal.sd = new double[Y.rows()];
  anal.doses = new double[Y.rows()];
  anal.prior = NULL;
  anal.isIncreasing = is_increasing;
  anal.alpha = 0.005; // alpha for analyses;
  anal.BMD_type = riskType;
  anal.BMR = bmrf;
  anal.samples = samples;
  anal.burnin = burnin;
  anal.tail_prob = tail_p;
  anal.suff_stat = Y.cols() == 3;

  for (int i = 0; i < Y.rows(); i++)
  {
    anal.Y[i] = Y(i, 0);
    anal.doses[i] = X(i, 0);
    if (Y.cols() == 3)
    { // sufficient statistics
      anal.n_group[i] = Y(i, 1);
      anal.sd[i] = Y(i, 2);
    }
  }

  estimate_ma_MCMC(&ma_anal, &anal, ma_result, &model_mcmc_info);

  List t1 = convert_mcmc_results(&model_mcmc_info);
  List t2 = convert_continuous_maresults_to_list(ma_result);
  List rV = List::create(Named("mcmc_runs") = t1, Named("ma_results") = t2);
  //////////////////////////////////////////////////////////
  // free up memory
  //////////////////////////////////////////////////////////
  for (int i = 0; i < ma_result->nmodels; i++)
  {
    del_continuous_model_result(ma_result->models[i]);
    del_mcmc_analysis(model_mcmc_info.analyses[i]);
  }
  delete[] ma_result->models;
  delete[] model_mcmc_info.analyses;
  delete[] ma_result->post_probs;
  delete[] ma_result->bmd_dist;
  delete ma_result;

  del_continuous_analysis(anal);
  del_continuousMA_analysis(ma_anal);
  return rV;
}

////////////////////////////////////////////////////////////////////////
// function: List run_ma_dichotomous()
// Purpose:  runs a model average based on the prior
//
// [[Rcpp::export(.run_ma_dichotomous)]]
List run_ma_dichotomous(Eigen::MatrixXd data, List priors, NumericVector models,
                        NumericVector model_p, bool is_MCMC,
                        NumericVector options1, IntegerVector options2)
{

  dichotomous_analysis Anal;
  Anal.BMD_type = (options2[2] == 1) ? eExtraRisk : eAddedRisk;
  Anal.BMR = options1[0];
  Anal.alpha = options1[1];
  Anal.Y = new double[data.rows()];
  Anal.n_group = new double[data.rows()];
  Anal.doses = new double[data.rows()];
  Anal.n = data.rows();
  Anal.samples = options2[2];
  Anal.burnin = options2[3];

  for (int i = 0; i < data.rows(); i++)
  {
    Anal.Y[i] = data(i, 1);
    Anal.n_group[i] = data(i, 2);
  }

  for (int i = 0; i < data.rows(); i++)
  {
    Anal.doses[i] = data(i, 0);
  }

  dichotomousMA_analysis ma_info;
  ma_info.nmodels = priors.size();
  ma_info.priors = new double *[priors.size()];
  ma_info.actual_parms = new int[priors.size()]; // actual number of parameters in the model
  ma_info.prior_cols = new int[priors.size()];
  ma_info.models = new int[priors.size()]; // given model
  ma_info.modelPriors = new double[priors.size()];

  for (int i = 0; i < priors.size(); i++)
  {

    Eigen::MatrixXd temp_cov = priors[i];

    ma_info.priors[i] = new double[temp_cov.rows() * temp_cov.cols()];
    ma_info.actual_parms[i] = temp_cov.rows();
    ma_info.prior_cols[i] = temp_cov.cols();
    ma_info.models[i] = models[i];
    ma_info.modelPriors[i] = model_p[i]; // prior over the model

    for (int m = 0; m < temp_cov.rows(); m++)
    {
      for (int n = 0; n < temp_cov.cols(); n++)
      {
        ma_info.priors[i][m + n * temp_cov.rows()] = temp_cov(m, n);
      }
    }
  }

  ma_MCMCfits model_mcmc_info;
  model_mcmc_info.nfits = ma_info.nmodels;
  model_mcmc_info.analyses = new bmd_analysis_MCMC *[ma_info.nmodels];
  dichotomousMA_result *ma_res = new_dichotomousMA_result(ma_info.nmodels, 300);
  for (int i = 0; i < ma_info.nmodels; i++)
  {
    // add a new result for each model result

    ma_res->models[i] = new_dichotomous_model_result(ma_info.models[i],
                                                     ma_info.actual_parms[i], 300);

    model_mcmc_info.analyses[i] = new_mcmc_analysis(ma_info.models[i],
                                                    ma_info.prior_cols[i],
                                                    Anal.samples);
  }

  List returnV;
  /////////
  if (is_MCMC)
  {
    estimate_ma_MCMC(&ma_info,
                     &Anal,
                     ma_res,
                     &model_mcmc_info);
    returnV = convert_dichotomous_maresults_to_list(ma_res);
    List t1 = convert_mcmc_results(&model_mcmc_info);
    returnV = List::create(Named("mcmc_runs") = t1, Named("ma_results") = returnV);
    // convert MCMC results to R Lists
    // List::create(Named("mcmc_runs") = t1 , Named("ma_results") = t2);
  }
  else
  {
    estimate_ma_laplace(&ma_info,
                        &Anal,
                        ma_res);
    // convert Laplace results to R lists
    returnV = convert_dichotomous_maresults_to_list(ma_res);
  }
  ///////////////////
  // to do add degree to the individual model
  for (int i = 0; i < priors.size(); i++)
  {
    delete[] ma_info.priors[i];
    del_mcmc_analysis(model_mcmc_info.analyses[i]);
  }

  delete[] model_mcmc_info.analyses;
  delete[] ma_info.priors;
  delete[] Anal.Y;
  delete[] Anal.n_group;
  delete[] Anal.doses;
  delete[] ma_info.actual_parms; // actual number of parameters in the model
  delete[] ma_info.prior_cols;
  delete[] ma_info.models;
  delete[] ma_info.modelPriors;
  delete_dichotomousMA_result(ma_res);

  return returnV;
}
