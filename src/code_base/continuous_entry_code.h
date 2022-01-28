/*
 * 
 * Copyright 2020  NIEHS <matt.wheeler@nih.gov>
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

#if defined R_COMPILATION || defined __cplusplus
  #include <cmath> 
#else
  #include <math.h>
#endif

#ifdef R_COMPILATION  
  #include <RcppEigen.h>
  #include <RcppGSL.h>
  using namespace Rcpp;
  using Rcpp::as;

#elif __cplusplus
  #include <Eigen/Dense>
#endif

#include "bmdStruct.h"

#if defined R_COMPILATION || defined __cplusplus
#include <gsl/gsl_randist.h>
#include <string>
#include <vector>
#include <limits>
//#include <math.h>
#include <stdio.h>
#include "bmds_entry.h"
#include "continuous_model_functions.h"


using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;

//const Map<MatrixXd> A(as<Map<MatrixXd>>(AA));

#include  <statmod.h>

#include <log_likelihoods.h>
#include <normal_likelihoods.h>
#include <normalModels.h>
#include <binomModels.h>
#include <IDPrior.h>

#include <bmd_calculate.h>

#include "normal_FUNL_NC.h"
#include "normal_HILL_NC.h"
#include "normal_POWER_NC.h"
#include "normal_POLYNOMIAL_NC.h"
#include "normal_EXP_NC.h"

#include "lognormal_HILL_NC.h"
#include "lognormal_POWER_NC.h"
#include "lognormal_POLYNOMIAL_NC.h"
#include "lognormal_EXP_NC.h"

#include "continuous_clean_aux.h"
#include "mcmc_analysis.h"

#endif

#ifndef _CONTINUOUS_ENTRY_CODE_H
#define _CONTINUOUS_ENTRY_CODE_H

#if defined R_COMPILATION || __cplusplus

template <class LL> 
Eigen::MatrixXd X_gradient_cont_norm(Eigen::MatrixXd theta, Eigen::MatrixXd Y, Eigen::MatrixXd D,bool SS,
                                     bool CV, int junk = 1){
  
  LL data_likelihood(Y,D,SS,CV,junk); 
  Eigen::MatrixXd rValue(Y.rows(),data_likelihood.nParms()) ;
  
  double *grad = new double[data_likelihood.nParms()+10]; 
  
  Eigen::MatrixXd md = D; 
  for (int i = 0; i < D.rows(); i++){
    md = D.row(i); 
    xgrad<LL>(theta, grad, &data_likelihood, md); 
    for (int j = 0; j < data_likelihood.nParms(); j++){
      rValue(i,j) =  grad[j]; 
    }
  }
  delete[] grad; 
  return rValue; 
  
}


template <class LL> 
Eigen::MatrixXd X_compute_mean_cont_norm(Eigen::MatrixXd theta, Eigen::MatrixXd Y, Eigen::MatrixXd D,bool SS,
                                         bool CV, int junk = 1){
  
  LL data_likelihood(Y,D,SS,CV,junk);
  Eigen::MatrixXd  md = D;
  
  Eigen::MatrixXd rValue = data_likelihood.mean(theta,md); 
  return rValue; 
}

template <class LL> 
Eigen::MatrixXd X_cov_cont_norm(Eigen::MatrixXd theta, Eigen::MatrixXd Y, Eigen::MatrixXd D,bool SS,
                                bool CV, int junk = 1){
  
  LL data_likelihood(Y,D,SS,CV,junk);
  
  Eigen::MatrixXd rValue = data_likelihood.variance(theta,D); //*Y.col(1).array(); // n*p*(1-p) 
  if (SS){
    rValue = (1/rValue.array())*Y.col(2).array(); // make up for the fact it is the variance
    return rValue.asDiagonal();                                         // for multiple observations
  }
  
  return rValue.asDiagonal().inverse(); 
  
}


template <class LL> 
Eigen::MatrixXd X_gradient_cont( Eigen::MatrixXd theta,Eigen::MatrixXd Y,
                                 Eigen::MatrixXd D,bool SS, int degree = 1){
  
  LL data_likelihood(Y,D,SS,degree); 
  Eigen::MatrixXd rValue(Y.rows(),data_likelihood.nParms()) ;
  
  double *grad = new double[data_likelihood.nParms()+10]; 
  
  Eigen::MatrixXd md = D; 
  for (int i = 0; i < D.rows(); i++){
    md = D.row(i); 
    xgrad<LL>(theta, grad, &data_likelihood, md); 
    for (int j = 0; j < data_likelihood.nParms(); j++){
      rValue(i,j) =  grad[j];   
    }
  }
  
  delete grad; 
  return rValue; 
  
}

template <class LL> 
Eigen::MatrixXd X_compute_mean_cont( Eigen::MatrixXd theta, Eigen::MatrixXd Y, Eigen::MatrixXd D,  
                                      bool SS, int degree = 1){
  
  LL data_likelihood(Y,D,SS,degree);
  Eigen::MatrixXd  md = D;
  Eigen::MatrixXd rValue = data_likelihood.mean(theta,md); 
  return rValue; 
}

template <class LL> 
Eigen::MatrixXd X_cov_cont( Eigen::MatrixXd theta,Eigen::MatrixXd Y, Eigen::MatrixXd D,  
                           bool SS, int degree = 1){
  
  LL data_likelihood(Y,D,SS,degree);
  Eigen::MatrixXd rValue = data_likelihood.variance(theta,D); //*Y.col(1).array(); // n*p*(1-p) 
 
  if (SS){
    rValue = (1/rValue.array())*Y.col(2).array(); // make up for the fact it is the variance
    return rValue.asDiagonal();                                         // for multiple observations
  }
  
  return rValue.asDiagonal().inverse(); 
  
}

bmd_analysis create_bmd_analysis_from_mcmc(unsigned int burnin, mcmcSamples s,double );
void transfer_mcmc_output(mcmcSamples a, bmd_analysis_MCMC *b); 

bool convertSStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                  Eigen::MatrixXd *SSTAT, Eigen::MatrixXd *SSTAT_LN,
                  Eigen::MatrixXd *UX);
                  
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeCol(Eigen::MatrixXd& matrix, unsigned int colToRemove);

bmd_analysis laplace_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                               Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                               bool is_increasing, 
                               double bmrf,   double bk_prob, 
                               double alpha, double step_size,
                               Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1,1), bool isFast = false);
                               
bmd_analysis laplace_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                            Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                            bool is_increasing, bool bConstVar,
                            double bmrf,   double bk_prob, 
                            double alpha, double step_size,
                            Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1,1),
                            int degree = 2, bool isFast = false);

                            
void transfer_continuous_model(bmd_analysis a, continuous_model_result *model); 

void bmd_range_find(continuousMA_result *res, 
					double *range);
					

					
/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a continuous Model Average (MA) for dichotomous
 *           data.  This is done using MCMC estimation. 
 * Input  : 
 *  @continuousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @continuous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @continuousMA_result   - Pointer to the MA result.  Inside the result are all of the 
 *                             individual fit information.  The fit information is returned
 *                             using a summary statistic format, and is computed using the
 *                             proper number of burnins. 
 *  @ma_MCMCfits            - Pointer to a structure that returns the individual MCMC results
 *   
 
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks, and may be problematic on 
 *                  server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data. 
 */					
void estimate_ma_MCMC(continuousMA_analysis    *MA ,
                           continuous_analysis *CA ,
                           continuousMA_result *res,
                           ma_MCMCfits         *ma);

/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a continuous Model Average (MA) for dichotomous
 *           data.  This is done using Laplace/MAP estimation. 
 * Input  : 
 *  @continuousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @continuous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @continuousMA_result   - Pointer to the MA result.  Inside the result are all of the 
 *                             individual fit information.  The fit information is returned
 *                             using a summary statistic format, and is computed using the
 *                             proper number of burnins. 
 *                             
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks, and may be problematic on 
 *                  server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data. 
 */						
void estimate_ma_laplace(continuousMA_analysis *MA,
                         continuous_analysis *CA ,
                         continuousMA_result *res);

/* Function: estimate_sm_laplace
 * Purpose:  This function performs a single model estimate for continuous data  
 * This is done using Laplace estimation. 
 * Input  : 
 *
 *  @continous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @continuous_model_result   - Pointer to a single model. The fit information is returned
 *                             using a summary statistic format.  All statistics are computed using
 *                             methodologies described in Wheeler et al (2020).
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */

void estimate_sm_laplace(continuous_analysis *CA ,
                         continuous_model_result *res, bool isFast = false);

/* Function: estimate_sm_mcmc
 * Purpose:  This function performs a single model estimate for continuous data  
 * This is done using MCMC estimation. 
 * Input  : 
 *
 *  @continous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @continuous_model_result   - Pointer to a single model. The fit information is returned
 *                             using a summary statistic format.  All statistics are computed using
 *                             methodologies described in Wheeler et al (2020).
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_sm_mcmc(continuous_analysis *CA,
                      continuous_model_result *res,
                      bmd_analysis_MCMC  *mcmc); 


/**************************************************
 * 
 * 
 ***************************************************/

void estimate_log_normal_aod(continuous_analysis *CA,
                            continuous_deviance *aod);

void estimate_normal_aod(continuous_analysis *CA,
                             continuous_deviance *aod);

void continuous_expectation( const continuous_analysis *CA, const continuous_model_result *MR,
                             continuous_expected_result *expected); 

#endif

#ifdef __cplusplus
extern "C" {
#endif
void estimate_sm_laplace_cont(struct continuous_analysis *CA,
                         struct continuous_model_result *res);
#ifdef __cplusplus
}
#endif

#endif
