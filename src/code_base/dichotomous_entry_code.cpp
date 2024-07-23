#include <cmath>
#ifdef R_COMPILATION
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
#include <gsl/gsl_randist.h>
// #include <math.h>
#include <stdio.h>

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "DichGammaBMD_NC.h"
#include "DichHillBMD_NC.h"
#include "DichLogLogisticBMD_NC.h"
#include "DichLogProbitBMD_NC.h"
#include "DichLogisticBMD_NC.h"
#include "DichMultistageBMD_NC.h"
#include "DichProbitBMD_NC.h"
#include "DichQlinearBMD_NC.h"
#include "DichWeibullBMD_NC.h"
#include "IDPrior.h"
#include "IDPriorMCMC.h"
#include "binomialTests.h"

#include "bmdStruct.h"
#include "bmds_entry.h"

#include "bmd_calculate.h"
#include "dichotomous_entry_code.h"
#include "mcmc_analysis.h"

using namespace Eigen;
void bmd_range_find(dichotomousMA_result *res, double *range) {
  // assume the minimum BMD for the MA is always 0
  range[0] = 0.0;
  double current_max = 0.0;
  for (int j = 10; j > 1; j--) {
    for (int i = 0; i < res->nmodels; i++) {
      int temp_idx = res->models[i]->dist_numE - j;

      // make sure we are not dealing with an infinite value
      // or not a number
      if (!isnan(res->models[i]->bmd_dist[temp_idx]) &&
          !isinf(res->models[i]->bmd_dist[temp_idx])) {
        if (res->models[i]->bmd_dist[temp_idx] > current_max) {
          current_max = res->models[i]->bmd_dist[temp_idx];
        }
      }
    }
  }
  // if we don't find a max then the upper limit is NAN
  range[1] = current_max == 0.0 ? std::numeric_limits<double>::quiet_NaN()
                                : current_max;
}

void rescale(Eigen::MatrixXd *parms, dich_model model, double max_dose) {

  Eigen::MatrixXd temp = *parms;
  Eigen::MatrixXd add = Eigen::MatrixXd::Zero(parms->rows(), 1);
  Eigen::MatrixXd scale = Eigen::MatrixXd::Ones(parms->rows(), 1);
  switch (model) {
  case dich_model::d_hill:
    add(2, 0) = temp(3, 0) * log(1 / max_dose);
    break;
  case dich_model::d_gamma:
    scale(2, 0) = 1 / max_dose;
    break;
  case dich_model::d_logistic:
    scale(1, 0) = 1 / max_dose;
    break;
  case dich_model::d_loglogistic:
    add(1, 0) = temp(2, 0) * log(1 / max_dose);
    break;
  case dich_model::d_logprobit:
    add(1, 0) = temp(2, 0) * log(1 / max_dose);
    break;
  case dich_model::d_multistage:
    for (int i = 1; i < parms->rows(); i++) {
      scale(i, 0) = pow(1 / max_dose, i);
    }
    break;
  case dich_model::d_probit:
    scale(1, 0) = 1 / max_dose;
    break;
  case dich_model::d_qlinear:
    scale(1, 0) = 1 / max_dose;
    break;
  case dich_model::d_weibull:
    scale(2, 0) = pow(1 / max_dose, temp(1, 0));
    break;
  }
  *parms = parms->array() * scale.array();
  *parms += add;
}

void rescale_var_matrix(Eigen::MatrixXd *var, Eigen::MatrixXd parms,
                        dich_model model, double max_dose) {
  Eigen::MatrixXd tvar = *var;
  Eigen::MatrixXd temp = parms;
  Eigen::MatrixXd add = Eigen::MatrixXd::Zero(parms.rows(), parms.rows());
  Eigen::MatrixXd scale = Eigen::MatrixXd::Identity(parms.rows(), parms.rows());

  switch (model) {
  case dich_model::d_hill:
    scale(2, 3) = log(1 / max_dose);
    break;
  case dich_model::d_gamma:
    scale(2, 2) = 1 / max_dose;
    break;
  case dich_model::d_logistic:
    scale(1, 1) = 1 / max_dose;
    break;
  case dich_model::d_loglogistic:
    scale(1, 2) = log(1 / max_dose);
    break;
  case dich_model::d_logprobit:
    scale(1, 2) = log(1 / max_dose);
    break;
  case dich_model::d_multistage:
    for (int i = 1; i < parms.rows(); i++) {
      scale(i, 0) = pow(1 / max_dose, i);
    }
    break;
  case dich_model::d_probit:
    scale(1, 1) = 1 / max_dose;
    break;
  case dich_model::d_qlinear:
    scale(1, 1) = 1 / max_dose;
    break;
  case dich_model::d_weibull:
    scale(2, 1) = log(1 / max_dose) * pow(1 / max_dose, temp(1, 0));
    break;
  }

  *var = scale * (tvar)*scale.transpose();
}
void rescale_dichotomous_model(mcmcSamples *v, dich_model model,
                               double max_dose) {
  // rescale the dichotomous to the origional values

  for (int i = 0; i < v->samples.cols(); i++) {
    Eigen::MatrixXd temp = v->samples.col(i);
    rescale(&temp, model, max_dose);
    v->samples.col(i) = temp;
    v->BMD(0, i) *= max_dose;
  }

  rescale_var_matrix(&v->map_cov, v->map_estimate, (dich_model)model, max_dose);
  rescale(&v->map_estimate, (dich_model)model, max_dose);
}

void transfer_dichotomous_model(bmd_analysis a,
                                dichotomous_model_result *model) {
  if (model) {
    model->nparms = a.COV.rows();
    model->max = a.MAP;
    model->bmd = a.MAP_BMD;
    for (int i = 0; i < model->dist_numE; i++) {
      double temp = double(i) / double(model->dist_numE);
      model->bmd_dist[i] = a.BMD_CDF.inv(temp); // BMD @ probability

      model->bmd_dist[model->dist_numE + i] = temp; // probability
    }

    for (int i = 0; i < model->nparms; i++) {
      model->parms[i] = a.MAP_ESTIMATE(i, 0);
      for (int j = 0; j < model->nparms; j++) {
        model->cov[i + j * model->nparms] = a.COV(i, j);
      }
    }
  }
}

void estimate_sm_mcmc(dichotomous_analysis *DA, dichotomous_model_result *res,
                      bmd_analysis_MCMC *mcmc, bool do_a_rescale) {

  ///////////////////////////////////
  Eigen::MatrixXd Y(DA->n, 2);
  Eigen::MatrixXd D(DA->n, 1);
  Eigen::MatrixXd prior(DA->parms, DA->prior_cols);
  for (int i = 0; i < DA->n; i++) {
    Y(i, 0) = DA->Y[i];
    Y(i, 1) = DA->n_group[i];
    D(i, 0) = DA->doses[i];
  }

  double max_dose = D.maxCoeff();
  D = (1 / max_dose) * D;

  for (int i = 0; i < DA->parms; i++) {
    for (int j = 0; j < DA->prior_cols; j++) {
      prior(i, j) = DA->prior[i + j * DA->parms];
    }
  }

  // copy the prior over.
  mcmcSamples a;
  std::vector<bool> fixedB;
  std::vector<double> fixedV;
  for (int i = 0; i < prior.rows(); i++) {
    fixedB.push_back(false);
    fixedV.push_back(0.0);
  }

  switch (DA->model) {
  case dich_model::d_hill:
    a = MCMC_bmd_analysis_DNC<dich_hillModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);

    break;
  case dich_model::d_gamma:
    a = MCMC_bmd_analysis_DNC<dich_gammaModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);

    break;
  case dich_model::d_logistic:
    a = MCMC_bmd_analysis_DNC<dich_logisticModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
    break;
  case dich_model::d_loglogistic:
    a = MCMC_bmd_analysis_DNC<dich_loglogisticModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
    break;
  case dich_model::d_logprobit:
    a = MCMC_bmd_analysis_DNC<dich_logProbitModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
    break;
  case dich_model::d_multistage:

    a = MCMC_bmd_analysis_DNC<dich_multistageNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
    break;
  case dich_model::d_probit:
    a = MCMC_bmd_analysis_DNC<dich_probitModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
    break;
  case dich_model::d_qlinear:
    a = MCMC_bmd_analysis_DNC<dich_qlinearModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
    break;
  case dich_model::d_weibull:
    a = MCMC_bmd_analysis_DNC<dich_weibullModelNC, IDPriorMCMC>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha, DA->samples);
  default:
    break;
  }

  if (do_a_rescale) {
  }

  bmd_analysis b;
  b = create_bmd_analysis_from_mcmc(DA->burnin, a, 1.0);
  mcmc->model = DA->model;
  mcmc->burnin = DA->burnin;
  mcmc->samples = DA->samples;
  mcmc->nparms = DA->parms;
  // rescale mcmc information

  rescale_dichotomous_model(&a, (dich_model)DA->model, max_dose);

  transfer_mcmc_output(a, mcmc);

  transfer_dichotomous_model(b, res);
  // rescale the BMD
  // has nothing to do with the do_a_rescale
  for (int i = 0; i < res->dist_numE; i++) {
    res->bmd_dist[i] *= max_dose;
  }
  res->bmd *= max_dose;
  res->model = DA->model;
  return;
}

void estimate_sm_laplace(dichotomous_analysis *DA,
                         dichotomous_model_result *res, bool do_a_rescale) {

  ///////////////////////////////////
  Eigen::MatrixXd Y(DA->n, 2);
  Eigen::MatrixXd D(DA->n, 1);
  Eigen::MatrixXd prior(DA->parms, DA->prior_cols);
  for (int i = 0; i < DA->n; i++) {
    Y(i, 0) = DA->Y[i];
    Y(i, 1) = DA->n_group[i];
    D(i, 0) = DA->doses[i];
  }

  double max_dose = D.maxCoeff();
  D = (1 / max_dose) * D;

  for (int i = 0; i < DA->parms; i++) {
    for (int j = 0; j < DA->prior_cols; j++) {
      prior(i, j) = DA->prior[i + j * DA->parms];
    }
  }
  // copy the prior over.
  bmd_analysis a;
  std::vector<bool> fixedB;
  std::vector<double> fixedV;
  for (int i = 0; i < DA->parms; i++) {
    fixedB.push_back(false);
    fixedV.push_back(0.0);
  }

  Eigen::MatrixXd Xd;
  Eigen::MatrixXd cv_t;
  Eigen::MatrixXd pr;

  switch ((dich_model)DA->model) {
  case dich_model::d_hill:
    a = bmd_analysis_DNC<dich_hillModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha * 0.5, 0.02);
    Xd = X_gradient<dich_hillModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_hillModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();

    break;
  case dich_model::d_gamma:

    a = bmd_analysis_DNC<dich_gammaModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha * 0.5, 0.02);
    Xd = X_gradient<dich_gammaModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_gammaModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();

    break;
  case dich_model::d_logistic:
    a = bmd_analysis_DNC<dich_logisticModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha * 0.5, 0.02);
    Xd = X_gradient<dich_logisticModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_logisticModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();
    break;
  case dich_model::d_loglogistic:
    a = bmd_analysis_DNC<dich_loglogisticModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha * 0.5, 0.02);
    Xd = X_gradient<dich_loglogisticModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_loglogisticModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();

    break;
  case dich_model::d_logprobit:
    a = bmd_analysis_DNC<dich_logProbitModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        DA->alpha * 0.5, 0.02);
    Xd = X_gradient<dich_logProbitModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_logProbitModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();

    break;
  case dich_model::d_multistage:
    a = bmd_analysis_DNC<dich_multistageNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        0.5 * DA->alpha, 0.02);

    Xd = X_gradient<dich_multistageNC>(a.MAP_ESTIMATE, Y, D, DA->degree);
    cv_t = X_cov<dich_multistageNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();

    break;
  case dich_model::d_probit:
    a = bmd_analysis_DNC<dich_probitModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        0.5 * DA->alpha, 0.02);
    Xd = X_gradient<dich_probitModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_probitModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();
    break;
  case dich_model::d_qlinear:
    a = bmd_analysis_DNC<dich_qlinearModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        0.5 * DA->alpha, 0.02);
    Xd = X_gradient<dich_qlinearModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_qlinearModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();
    break;
  case dich_model::d_weibull:

    a = bmd_analysis_DNC<dich_weibullModelNC, IDPrior>(
        Y, D, prior, fixedB, fixedV, DA->degree, DA->BMR, DA->BMD_type,
        0.5 * DA->alpha, 0.02);
    Xd = X_gradient<dich_weibullModelNC>(a.MAP_ESTIMATE, Y, D);
    cv_t = X_cov<dich_weibullModelNC>(a.MAP_ESTIMATE, Y, D);
    pr = X_logPrior<IDPrior>(a.MAP_ESTIMATE, prior);
    pr = Xd.transpose() * cv_t * Xd + pr;
    Xd = Xd * pr.inverse() * Xd.transpose() * cv_t;
    res->model_df = Xd.diagonal().array().sum();
  default:
    break;
  }

  if (do_a_rescale) {
    rescale_var_matrix(&a.COV, a.MAP_ESTIMATE, (dich_model)DA->model, max_dose);

    rescale(&a.MAP_ESTIMATE, (dich_model)DA->model, max_dose);
  }

  transfer_dichotomous_model(a, res);
  res->bmd *= max_dose;
  // rescale the BMD
  for (int i = 0; i < res->dist_numE; i++) {
    res->bmd_dist[i] *= max_dose;
  }

  res->model = DA->model;
  return;
}

void estimate_ma_MCMC(dichotomousMA_analysis *MA, dichotomous_analysis *DA,
                      dichotomousMA_result *res, ma_MCMCfits *ma) {

#pragma omp parallel for
  for (int i = 0; i < MA->nmodels; i++) {
    dichotomous_analysis temp = *DA; // copy over the initial stuff
    temp.prior = MA->priors[i];
    temp.parms = MA->actual_parms[i];
    temp.prior_cols = MA->prior_cols[i];
    temp.model = MA->models[i];
    if (MA->models[i] == dich_model::d_multistage) {
      temp.degree = temp.parms - 1;
    } else {
      temp.degree = 0;
    }
    // fit the individual model
    estimate_sm_mcmc(&temp, res->models[i], ma->analyses[i], false);
  }
  double *post_probs = new double[MA->nmodels];
  double temp = 0.0;
  double max_prob = -1.0 * std::numeric_limits<double>::infinity();
  for (int i = 0; i < MA->nmodels; i++) {
    Eigen::Map<MatrixXd> transfer_mat(
        res->models[i]->cov, res->models[i]->nparms, res->models[i]->nparms);
    Eigen::MatrixXd cov = transfer_mat;
    temp = res->models[i]->nparms / 2 * log(2 * M_PI) - res->models[i]->max +
           0.5 * log(max(0.0, cov.determinant()));
    if (cov.determinant() < 0 || !std::isfinite(res->models[i]->bmd)) {
      temp = -1.0 * std::numeric_limits<double>::infinity();
    }
    if (isfinite(temp)) {
      max_prob = temp > max_prob ? temp : max_prob;
    } else {
      temp = -1.0 * std::numeric_limits<double>::infinity();
    }
    post_probs[i] = temp;
  }
  double max_dose = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < DA->n; i++) {
    if (DA->doses[i] > max_dose) {
      max_dose = DA->doses[i];
    }
  }

  // now that the Bayes Factor's have been approximated, we can rescale
  for (int i = 0; i < MA->nmodels; i++) {
    Eigen::Map<MatrixXd> trans_cov(res->models[i]->cov, res->models[i]->nparms,
                                   res->models[i]->nparms);
    Eigen::MatrixXd cov = trans_cov;
    Eigen::Map<MatrixXd> trans_par(res->models[i]->parms,
                                   res->models[i]->nparms, 1);
    Eigen::MatrixXd par = trans_par;

    rescale(&par, (dich_model)res->models[i]->model, max_dose);
    rescale_var_matrix(&cov, par, (dich_model)res->models[i]->model, max_dose);

    for (int n = 0; n < res->models[i]->nparms; n++) {
      res->models[i]->parms[n] = par(n, 0);
      for (int m = 0; m < res->models[i]->nparms; m++) {
        res->models[i]->cov[n + m * res->models[i]->nparms] = cov(n, m);
      }
    }
  }

  double norm_sum = 0.0;

  for (int i = 0; i < MA->nmodels; i++) {
    post_probs[i] = post_probs[i] - max_prob + log(MA->modelPriors[i]);
    norm_sum += exp(post_probs[i]);
    post_probs[i] = exp(post_probs[i]);
  }

  for (int j = 0; j < MA->nmodels; j++) {
    post_probs[j] = post_probs[j] / norm_sum;
    // cout << post_probs[j] << endl;

    for (int i = 0; i < res->models[j]->dist_numE; i++) {

      if (isnan(res->models[j]->bmd_dist[i])) {
        post_probs[j] =
            0; // if the cdf has nan in it then it needs a 0 posterior
      }
    }
  }

  norm_sum = 0.0;
  for (int i = 0; i < MA->nmodels; i++) {
    norm_sum += post_probs[i];
  }

  for (int i = 0; i < MA->nmodels; i++) {
    post_probs[i] = post_probs[i] / norm_sum;
    res->post_probs[i] = post_probs[i];
    res->models[i]->model = MA->models[i];
  }

  double range[2];
  std::vector<bmd_cdf> model_cdfs(MA->nmodels);

  for (int i = 0; i < MA->nmodels; i++) {
    std::vector<double> bmd;
    ;
    std::vector<double> prc;
    for (int m = 0; m < res->models[i]->dist_numE; m++) {
      if (!isinf(res->models[i]->bmd_dist[m]) &&
          !isnan(res->models[i]->bmd_dist[m])) {
        // deal with numerical problems in the tails of the distribution
        if (m > 0) {
          if (res->models[i]->bmd_dist[m] <= res->models[i]->bmd_dist[m - 1]) {
            res->models[i]->bmd_dist[m] =
                res->models[i]->bmd_dist[m - 1] + 1e-8;
          }
        }

        bmd.push_back(res->models[i]->bmd_dist[m]);
        prc.push_back(res->models[i]->bmd_dist[m + res->models[i]->dist_numE]);
      }
    }
    if (prc.size() > 5) {
      model_cdfs[i] = bmd_cdf(prc, bmd);
    }
  }

  bmd_range_find(res, range);
  double range_bmd = range[1] - range[0];

  for (int i = 0; i < res->dist_numE; i++) {
    double cbmd = double(i) / double(res->dist_numE) * range_bmd;
    double prob = 0.0;

    for (int j = 0; j < MA->nmodels; j++) {
      prob += isnan(model_cdfs[j].P(cbmd))
                  ? 0.0
                  : model_cdfs[j].P(cbmd) * post_probs[j];
    }
    // cout << prob << endl;
    res->bmd_dist[i] = cbmd;
    res->bmd_dist[i + res->dist_numE] = prob;
  }

  delete[] post_probs;
  return;
}

void estimate_ma_laplace(dichotomousMA_analysis *MA, dichotomous_analysis *DA,
                         dichotomousMA_result *res) {

#pragma omp parallel for
  for (int i = 0; i < MA->nmodels; i++) {
    dichotomous_analysis temp = *DA; // copy over the initial stuff
    temp.prior = MA->priors[i];
    temp.parms = MA->actual_parms[i];
    temp.prior_cols = MA->prior_cols[i];
    temp.model = MA->models[i];
    if (MA->models[i] == dich_model::d_multistage) {
      temp.degree = temp.parms - 1;
    } else {
      temp.degree = 0;
    }
    // fit the individual model
    estimate_sm_laplace(&temp, res->models[i], false);
  }

  double *post_probs = new double[MA->nmodels];
  double temp = 0.0;
  double max_prob = -1.0 * std::numeric_limits<double>::infinity();
  for (int i = 0; i < MA->nmodels; i++) {
    Eigen::Map<MatrixXd> transfer_mat(
        res->models[i]->cov, res->models[i]->nparms, res->models[i]->nparms);
    Eigen::MatrixXd cov = transfer_mat;
    temp = res->models[i]->nparms / 2 * log(2 * M_PI) - res->models[i]->max +
           0.5 * log(max(0.0, cov.determinant()));
    if (cov.determinant() < 0 || !std::isfinite(res->models[i]->bmd)) {
      temp = -1.0 * std::numeric_limits<double>::infinity();
    }
    if (isfinite(temp)) {
      max_prob = temp > max_prob ? temp : max_prob;
    } else {
      temp = -1.0 * std::numeric_limits<double>::infinity();
    }

    post_probs[i] = temp;
  }

  double max_dose = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < DA->n; i++) {
    if (DA->doses[i] > max_dose) {
      max_dose = DA->doses[i];
    }
  }
  // now that the Bayes Factor's have been approximated, we can rescale
  for (int i = 0; i < MA->nmodels; i++) {
    Eigen::Map<MatrixXd> trans_cov(res->models[i]->cov, res->models[i]->nparms,
                                   res->models[i]->nparms);
    Eigen::MatrixXd cov = trans_cov;
    Eigen::Map<MatrixXd> trans_par(res->models[i]->parms,
                                   res->models[i]->nparms, 1);
    Eigen::MatrixXd par = trans_par;

    rescale(&par, (dich_model)res->models[i]->model, max_dose);
    rescale_var_matrix(&cov, par, (dich_model)res->models[i]->model, max_dose);

    for (int n = 0; n < res->models[i]->nparms; n++) {
      res->models[i]->parms[n] = par(n, 0);
      for (int m = 0; m < res->models[i]->nparms; m++) {
        res->models[i]->cov[n + m * res->models[i]->nparms] = cov(n, m);
      }
    }
  }

  double norm_sum = 0.0;

  for (int i = 0; i < MA->nmodels; i++) {
    post_probs[i] = post_probs[i] - max_prob + log(MA->modelPriors[i]);
    norm_sum += exp(post_probs[i]);
    post_probs[i] = exp(post_probs[i]);
  }

  for (int j = 0; j < MA->nmodels; j++) {
    post_probs[j] = post_probs[j] / norm_sum;
    for (int i = 0; i < res->models[j]->dist_numE; i++) {
      if (isnan(res->models[j]->bmd_dist[i])) {
        post_probs[j] =
            0; // if the cdf has nan in it then it needs a 0 posterior
      }
    }
  }

  norm_sum = 0.0;
  for (int i = 0; i < MA->nmodels; i++) {
    norm_sum += post_probs[i];
  }

  for (int i = 0; i < MA->nmodels; i++) {
    post_probs[i] = post_probs[i] / norm_sum;
    res->post_probs[i] = post_probs[i];
    res->models[i]->model = MA->models[i];
  }

  double range[2];

  std::vector<bmd_cdf> model_cdfs(MA->nmodels);
  for (int i = 0; i < MA->nmodels; i++) {
    std::vector<double> bmd;
    ;
    std::vector<double> prc;
    for (int m = 0; m < res->models[i]->dist_numE; m++) {
      if (!isinf(res->models[i]->bmd_dist[m]) &&
          !isnan(res->models[i]->bmd_dist[m])) {
        // deal with numerical problems in the tails
        if (m > 0) {
          if (res->models[i]->bmd_dist[m] <= res->models[i]->bmd_dist[m - 1]) {
            res->models[i]->bmd_dist[m] =
                res->models[i]->bmd_dist[m - 1] + 1e-8;
          }
        }

        bmd.push_back(res->models[i]->bmd_dist[m]);
        prc.push_back(res->models[i]->bmd_dist[m + res->models[i]->dist_numE]);
      }
    }
    if (prc.size() > 5) {
      model_cdfs[i] = bmd_cdf(prc, bmd);
    }
  }

  bmd_range_find(res, range);

  double range_bmd = range[1] - range[0];

  for (int i = 0; i < res->dist_numE; i++) {
    double cbmd = double(i) / double(res->dist_numE) * range_bmd;
    double prob = 0.0;

    for (int j = 0; j < MA->nmodels; j++) {
      prob += isnan(model_cdfs[j].P(cbmd))
                  ? 0.0
                  : model_cdfs[j].P(cbmd) * post_probs[j];
    }
    res->bmd_dist[i] = cbmd;
    res->bmd_dist[i + res->dist_numE] = prob;
  }
  delete[] post_probs;
  return;
}

void compute_dichotomous_pearson_GOF(dichotomous_PGOF_data *data,
                                     dichotomous_PGOF_result *res) {
  Eigen::MatrixXd Y(data->n, 2);
  Eigen::MatrixXd D(data->n, 1);
  Eigen::MatrixXd parms(data->parms, 1);
  Eigen::MatrixXd mean_p;
  Eigen::MatrixXd mean_d;

  for (int i = 0; i < data->parms; i++) {
    parms(i, 0) = data->est_parms[i];
  }

  for (int i = 0; i < data->n; i++) {
    Y(i, 0) = data->Y[i];
    Y(i, 1) = data->n_group[i];
    D(i, 0) = data->doses[i];
  }
  int degree = 1;
  switch (data->model) {
  case dich_model::d_hill:
    mean_d = X_compute_mean<dich_hillModelNC>(Y, D, parms);
    break;
  case dich_model::d_gamma:
    mean_d = X_compute_mean<dich_gammaModelNC>(Y, D, parms);
    break;
  case dich_model::d_logistic:
    mean_d = X_compute_mean<dich_logisticModelNC>(Y, D, parms);
    break;
  case dich_model::d_loglogistic:
    mean_d = X_compute_mean<dich_loglogisticModelNC>(Y, D, parms);
    break;
  case dich_model::d_logprobit:
    mean_d = X_compute_mean<dich_logProbitModelNC>(Y, D, parms);
    break;
  case dich_model::d_multistage:
    degree = parms.rows() - 1;
    mean_d = X_compute_mean<dich_multistageNC>(Y, D, parms, degree);
    break;
  case dich_model::d_probit:
    mean_d = X_compute_mean<dich_probitModelNC>(Y, D, parms);
    break;
  case dich_model::d_qlinear:
    mean_d = X_compute_mean<dich_qlinearModelNC>(Y, D, parms);
    break;
  case dich_model::d_weibull:
    mean_d = X_compute_mean<dich_weibullModelNC>(Y, D, parms);
  default:
    break;
  }

  Eigen::MatrixXd expected = Y.col(1).array() * mean_d.array();

  Eigen::MatrixXd residual = Y.col(0) - expected;
  residual = residual.array() / sqrt(expected.array());
  Eigen::MatrixXd sqresid = residual.array() * residual.array();
  Eigen::MatrixXd resultsTable(Y.rows(), 5);
  resultsTable << Y.col(0), Y.col(1), expected, residual, sqresid;
  for (int i = 0; i < data->n; i++) {
    res->expected[i] = expected(i, 0);
    res->residual[i] = residual(i, 0);
  }
  res->n = data->n; // total number of observations obs
  res->test_statistic = sqresid.array().sum();

  if (data->n - data->model_df > 0.0) {
    res->p_value =
        1.0 - gsl_cdf_chisq_P(sqresid.array().sum(), data->n - data->model_df);
  } else {
    res->p_value = 1.0;
  }
  res->df = data->n - data->model_df;
}

void estimate_sm_laplace_dicho(dichotomous_analysis *DA,
                               dichotomous_model_result *res,
                               bool do_a_rescale) {

  estimate_sm_laplace(DA, res, do_a_rescale);
}

Eigen::MatrixXd A1_startingValues(Eigen::MatrixXd X, Eigen::MatrixXd Y) {

  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  std::vector<double> udoses = vec; // this should be the unique dose group

  Eigen::MatrixXd meanX = Eigen::MatrixXd::Zero(Y.rows(), udoses.size());

  for (int i = 0; i < meanX.rows(); i++) {
    for (size_t j = 0; j < udoses.size(); j++) {
      meanX(i, j) = udoses[j] == X(i, 0) ? Y(i, 1) : 0.0;
    }
  }
  // cout << meanX << endl << ":" << endl;
  Eigen::MatrixXd temp = meanX.transpose() * meanX;
  temp = temp.inverse() * meanX.transpose() * Y.col(0);
  // logit transform
  temp = log(temp.array() / (1 - temp.array()));
  // check for P = 0 or P =1 i.e. Infinity
  for (int i = 0; i < temp.rows(); i++) {
    if (isinf(temp(i, 0))) {
      if (temp(i, 0) < 0) {
        temp(i, 0) = -17;
      } else {
        temp(i, 0) = 17;
      }
    }
  }

  return temp;
}

Eigen::MatrixXd A2_startingValues(Eigen::MatrixXd X, Eigen::MatrixXd Y) {

  Eigen::MatrixXd meanX = Eigen::MatrixXd::Zero(Y.rows(), 1).array();

  for (int i = 0; i < meanX.rows(); i++) {
    meanX(i, 0) = Y(i, 1);
  }
  Eigen::MatrixXd temp = meanX.transpose() * meanX;
  temp = temp.inverse() * meanX.transpose() * Y.col(0);
  // logit transform
  temp = log(temp.array() / (1 - temp.array()));
  // check for P = 0 or P =1 i.e. Infinity
  for (int i = 0; i < temp.rows(); i++) {
    if (isinf(temp(i, 0))) {
      if (temp(i, 0) < 0) {
        temp(i, 0) = -17;
      } else {
        temp(i, 0) = 17;
      }
    }
  }
  return temp;
}

void deviance_dichotomous(dichotomous_analysis *DA, dichotomous_aod *AOD) {

  ///////////////////////////////////
  Eigen::MatrixXd Y(DA->n, 2);
  Eigen::MatrixXd D(DA->n, 1);
  for (int i = 0; i < DA->n; i++) {
    Y(i, 0) = DA->Y[i];
    Y(i, 1) = DA->n_group[i];
    D(i, 0) = DA->doses[i];
  }

  binomialLLTESTA1 like_A1(Y, D);
  binomialLLTESTA2 like_A2(Y, D);
  Eigen::MatrixXd P1(like_A1.nParms(), 5);
  for (int i = 0; i < P1.rows(); i++) {
    P1.row(i) << 0, 0, 1, -30, 30;
  }
  Eigen::MatrixXd P2(like_A2.nParms(), 5);
  for (int i = 0; i < P2.rows(); i++) {
    P2.row(i) << 0, 0, 1, -30, 30;
  }
  std::vector<double> fix1(like_A1.nParms());
  for (int i = 0; i < like_A1.nParms(); i++) {
    fix1[i] = 0.0;
  }
  std::vector<bool> isfix1(like_A1.nParms());
  for (int i = 0; i < like_A1.nParms(); i++) {
    isfix1[i] = false;
  }

  IDcontinuousPrior P1Init(P1);
  statModel<binomialLLTESTA1, IDcontinuousPrior> a1Model(like_A1, P1Init,
                                                         isfix1, fix1);
  A1_startingValues(D, Y);
  optimizationResult a1Result = findMAP<binomialLLTESTA1, IDcontinuousPrior>(
      &a1Model, A1_startingValues(D, Y), OPTIM_NO_FLAGS);

  AOD->A1 = a1Result.functionV;
  AOD->N1 = a1Result.max_parms.rows();

  std::vector<double> fix2(like_A2.nParms());
  for (int i = 0; i < like_A2.nParms(); i++) {
    fix1[i] = 0.0;
  }
  std::vector<bool> isfix2(like_A2.nParms());
  for (int i = 0; i < like_A2.nParms(); i++) {
    isfix1[i] = false;
  }

  IDcontinuousPrior P2Init(P2);
  statModel<binomialLLTESTA2, IDcontinuousPrior> a2Model(like_A2, P2Init,
                                                         isfix2, fix2);
  optimizationResult a2Result = findMAP<binomialLLTESTA2, IDcontinuousPrior>(
      &a2Model, A2_startingValues(D, Y), OPTIM_NO_FLAGS);
  AOD->A2 = a2Result.functionV;
  AOD->N2 = a2Result.max_parms.rows();
}
