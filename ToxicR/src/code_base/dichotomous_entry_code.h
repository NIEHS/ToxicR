#include "bmdStruct.h"
//#include <math.h>
#if defined R_COMPILATION || defined __cplusplus
  #include <cmath>
#else
  #include <math.h>
#endif
#include <stdio.h>


#if defined R_COMPILATION || defined __cplusplus
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

#include "DichHillBMD_NC.h"
#include "DichMultistageBMD_NC.h"
#include "DichLogLogisticBMD_NC.h"
#include "DichLogProbitBMD_NC.h"
#include "DichWeibullBMD_NC.h"
#include "DichGammaBMD_NC.h"
#include "DichQlinearBMD_NC.h"
#include "DichLogisticBMD_NC.h"
#include "DichProbitBMD_NC.h"
#include "IDPrior.h"
#include "gradient.h"

#include "bmds_entry.h"
//#include "bmdStruct.h"
#include "continuous_entry_code.h"
#endif

#ifndef _DICHOTOMOUS_ENTRY_CODE_H
#define _DICHOTOMOUS_ENTRY_CODE_H

#if (defined R_COMPILATION || defined __cplusplus) 

template <class LL> 
Eigen::MatrixXd X_gradient( Eigen::MatrixXd theta,Eigen::MatrixXd Y,
                                      Eigen::MatrixXd D,int degree = 1){
        
        LL data_likelihood(Y,D,degree); 
        Eigen::MatrixXd rValue(Y.rows(),data_likelihood.nParms()) ;
       
        double grad[data_likelihood.nParms()]; 
   
        Eigen::MatrixXd md; 
        for (int i = 0; i < D.rows(); i++){
                md = data_likelihood.convertDataMatrix(D.row(i)); 
                
                xgrad<LL>(theta, grad, &data_likelihood, md); 
                for (int j = 0; j < data_likelihood.nParms(); j++){
                        rValue(i,j) =  grad[j] *Y(i,1); //n*p'  
                }
        }
        
        return rValue; 
        
}

template <class LL> 
Eigen::MatrixXd X_compute_mean( Eigen::MatrixXd Y, Eigen::MatrixXd D,  
                                Eigen::MatrixXd parms, int degree = 1){
     
     
     LL data_likelihood(Y,D,degree);
     Eigen::MatrixXd  md = data_likelihood.convertDataMatrix(D);
     Eigen::MatrixXd rValue = data_likelihood.mean(parms,md); 
     return rValue; 
     
}
template <class LL> 
Eigen::MatrixXd X_cov( Eigen::MatrixXd theta,Eigen::MatrixXd Y,
                            Eigen::MatrixXd D){
        
        LL data_likelihood(Y,D,1); 
        
        Eigen::MatrixXd rValue = data_likelihood.variance(theta); //*Y.col(1).array(); // n*p*(1-p) 
        rValue = rValue.array()*Y.col(1).array(); 
        return rValue.asDiagonal().inverse(); 
        
}

void compute_dichotomous_pearson_GOF(dichotomous_PGOF_data *data,dichotomous_PGOF_result *res);

/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a dichotomous Model Average (MA) for dichotomous
 *           data.  This is done using MCMC estimation. 
 * Input  : 
 *  @dichotomousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @dichotomousMA_result   - Pointer to the MA result.  Inside the result are all of the 
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
void estimate_ma_MCMC(dichotomousMA_analysis *MA,
                      dichotomous_analysis   *DA,
                      dichotomousMA_result   *res,
                      ma_MCMCfits            *ma);

/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a dichotomous Model Average (MA) for dichotomous
 *           data.  This is done using Laplace estimation. 
 * Input  : 
 *  @dichotomousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @dichotomousMA_result   - Pointer to the MA result.  Inside the result are all of the 
 *                             individual fit information.  The fit information is returned
 *                             using a summary statistic format.  This is computed using the 
 *                             methodologies described in Wheeler et al ()2020).
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_ma_laplace(dichotomousMA_analysis *MA,
                         dichotomous_analysis *DA ,
                         dichotomousMA_result *res);

/* Function: estimate_sm_laplace
 * Purpose:  This function performs a single model estimate for dichotomous data  This is done using Laplace estimation. 
 * Input  : 
 *
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @dichotomous_model_result   - Pointer to a single model. The fit information is returned
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
void estimate_sm_laplace(dichotomous_analysis *DA ,
                         dichotomous_model_result *res,
                         bool do_a_rescale = true);



/* Function: estimate_sm_mcmc
 * Purpose:  This function performs a single model estimate for dichotomous data  
 * This is done using MCMC estimation. 
 * Input  : 
 *
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @dichotomous_model_result   - Pointer to a single model. The fit information is returned
 *                             using a summary statistic format.  All statistics are computed using
 *                             methodologies described in Wheeler et al (2020).
 *  
 *  @bmd_analysis_MCMC          - All of the posterior samples from the MCMC analysis. 
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_sm_mcmc(dichotomous_analysis *DA, 
                      dichotomous_model_result *res,
                      bmd_analysis_MCMC *mcmc,
                      bool do_a_rescale = true);

void deviance_dichotomous(dichotomous_analysis *DA,
                              dichotomous_aod *AOD); 

void estimate_normal_variance(continuous_analysis *CA,
                              double *v_c, double *v_nc, double *v_pow);

#endif

//c entry
#ifdef __cplusplus
extern "C" {
#endif
  void estimate_sm_laplace_dicho(struct dichotomous_analysis *DA,
                                 struct dichotomous_model_result *res,
                                 bool do_a_rescale);
#ifdef __cplusplus
}
#endif

#endif
