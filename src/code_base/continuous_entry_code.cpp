/*
 * continuous_entry_code.cpp
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


 /* 
// a column order matrix px5 where p is 
// the number of parametersd
*/

#include "continuous_entry_code.h"
#include "analysis_of_deviance.h"
#include "mcmc_analysis.h"
#include "bmd_calculate.h"
#include "IDPriorMCMC.h"
 
#include <chrono>
#include <algorithm>
#include <vector>
#include <limits>
#include <set>
#include <numeric>

using namespace std; 

Eigen::MatrixXd quadraticRegression(Eigen::MatrixXd Y_N, Eigen::MatrixXd X){
  
  Eigen::MatrixXd mX = Eigen::MatrixXd::Zero(Y_N.rows(), 3); 
  Eigen::MatrixXd W  = Eigen::MatrixXd::Zero(Y_N.rows(), Y_N.rows());
  for (int i = 0; i < mX.rows(); i++)
  {  
    W(i, i) = Y_N.cols() == 3? pow(1/Y_N(i,1),2) * Y_N(i,2) : 1;
    for (int j = 0; j < 3; j++) {
      switch(j){
      case 2:
        mX(i, j) = X(i,0)*X(i,0);
        break; 
      case 1:
        mX(i, j) = X(i,0); 
        break; 
      default: 
        mX(i, j) = 1 ;
      break;  
      } 
    }
  }

  Eigen::MatrixXd betas = mX.transpose()*W*mX;
  betas = betas.inverse()*mX.transpose()*W*Y_N.col(0);
  return betas;
}

Eigen::MatrixXd powerSearchRegression(Eigen::MatrixXd Y_N, Eigen::MatrixXd X){
  
  double min = 0; 
  double sum ; 
  
  Eigen::MatrixXd mX = Eigen::MatrixXd::Zero(Y_N.rows(), 2); 
  Eigen::MatrixXd W  = Eigen::MatrixXd::Zero(Y_N.rows(), Y_N.rows());
  Eigen::MatrixXd E  = mX; 
  Eigen::MatrixXd betas;
  Eigen::MatrixXd rbetas(3,1);
  
  for (double pows = 1.0; pows < 17; pows += 0.5){
 
    for (int i = 0; i < mX.rows(); i++)
    {  
      W(i, i) = Y_N.cols() == 3? pow(1/Y_N(i,1),2) * Y_N(i,2) : 1;
      for (int j = 0; j < 2; j++) {
        switch(j){
        case 1:
          mX(i, j) = pow(X(i,0),pows);
          break; 
        default: 
          mX(i, j) = 1 ;
          break;  
        } 
      }
    }
    
    betas = mX.transpose()*W*mX;
    betas = betas.inverse()*mX.transpose()*W*Y_N.col(0);
    
    E = (Y_N.col(0) - mX * betas);
    E = E.array() * E.array(); 
    E = W*E; 
    sum = E.array().sum();
    
    if (pows == 1.0 || sum < min){
      rbetas(0,0) = betas(0,0);  
      rbetas(1,0) = betas(1,0); 
      rbetas(2,0) = pows; 
      min    =  sum; 
    }

  }
  return rbetas;
}
////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd init_funl_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
  prior(0,1) = betas(0,0); 
  double max_d = vec[vec.size()-1]; 
  double max_r = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d);
  prior(1,1)   = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d - prior(0,1))/max_d; 
  prior(2,1)   = max_r; 
  prior(3,1)   = 0.5;
  prior(4,1)   = 1;
  prior(5,1)   = 0.75;
  prior(6,1)   = 1;
  
  for (int i = 0; i < 7; i++){
    if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3); 
    if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
  }
  
  
  return prior; 
}
///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd init_hill_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  double minDose = X.minCoeff(); 
  double maxDose = X.maxCoeff(); 
  double init = 0; 
  int nmin = 0, nmax = 0;  
  
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++; 
      init += Y_N(i,0); 
    }
  }
  init *= init/double(nmin); 
  
  Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
  prior(0,1) = init; 
  init = 0;
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==maxDose){
      nmax++; 
      init += Y_N(i,0); 
    }
  }
  init *= init/double(nmin); 
  
  prior(1,1)   =  (init - prior(0,1))/(maxDose-minDose); 
  prior(2,1)   = 0;//0.0001*maxDose; 
  prior(3,1)   = 10;
  
  if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
  if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
  if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
  if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  //cerr << prior << endl; 
  return prior; 
}


Eigen::MatrixXd init_pow_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd betas = powerSearchRegression(Y_N,  X);
  prior(0,1)   = betas(0,0); 
  prior(1,1)   = betas(1,0);  
  prior(2,1)   = betas(2,0);
  
  if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
  if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
  if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
  if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  
  if (prior(2,1) < prior(1,3)) prior(2,1) = prior(1,3); 
  if (prior(2,1) > prior(1,4)) prior(2,1) = prior(1,4);
  
  return prior; 
}

Eigen::MatrixXd init_hill_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  Y_LN.col(0) = exp(Y_LN.col(0).array());
  if (Y_LN.cols() ==3 ){
    Y_LN.col(1) = exp(Y_LN.col(1).array());
  }
  return init_hill_nor(Y_LN,  X, prior); 
  
}


Eigen::MatrixXd init_exp_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
  prior(0,1) = betas(0,0); 
  double max_d = vec[vec.size()-1]; 
  double max_r = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d);
  prior(2,1)   = log(0.001);  
  double temp = max_r/prior(0,1);
  
  temp =  -(temp-exp(prior(2,1)))/(exp(prior(2,1))-1.0);
  
  prior(1,1) = 0.05; 
  prior(3,1)   = 2.5; 
  
  if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
  if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
  if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
  if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  
  
  return prior; 
}

Eigen::MatrixXd init_exp_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
 //right here
  Y_LN.col(0) = exp(Y_LN.col(0).array());
  if (Y_LN.cols() ==3 ){
      Y_LN.col(1) = exp(Y_LN.col(1).array());
  }
  return init_exp_nor(Y_LN, X, prior); 
}

Eigen::MatrixXd init_poly(Eigen::MatrixXd Y, Eigen::MatrixXd tX, 
                          Eigen::MatrixXd prior, int deg = 2){

  Eigen::MatrixXd X = Eigen::MatrixXd::Ones(tX.rows(),deg+1);
  Eigen::MatrixXd W = Eigen::MatrixXd::Identity(tX.rows(),tX.rows());
 
  for (int i = 0; i < X.rows(); i++){
    if (Y.cols()>1){
      W(i,i) = Y(i,2)/Y(i,1)*Y(i,1); // Weights: \sigma^2/N
    }
    for (int j = 1; j < X.cols(); j++){
      X(i,j) = pow(tX(i,0),j);  
    }
  }
  Eigen::MatrixXd B = Eigen::MatrixXd::Ones(deg+1,1);
  B = X.transpose()*W*X;
  B = B.inverse()*X.transpose()*W*Y.col(0); 
  for(int i = 0; i < B.rows(); i++){
    if ( B(i,0) < prior(i,3) ){
      prior(i,1) = prior(i,3); 
    } else if (B(i,0) > prior(i,4)){
      prior(i,1) = prior(i,4); 
    }else{
      prior(i,1) = B(i,0);
    }
  } 
  
  return prior; 
}

/*initialize_mle
 * This function is for MLE optimization it takes the data/model type and then tries 
 * to start the initializer on reasonable initial values. These values will then be fed
 * to the optimizer.
 * OUTPUT: new prior vector with the new initial values put in the mean column. 
 */
Eigen::MatrixXd initialize_model(Eigen::MatrixXd Y_N, Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, 
                                 Eigen::MatrixXd prior, distribution data_dist, cont_model model){
  
  Eigen::MatrixXd retVal = prior; 
  int deg = 0; 
  switch (model){
  case cont_model::funl:
    retVal = init_funl_nor(Y_N, X, prior);
    break; 
  case cont_model::hill:
    retVal = distribution::log_normal==data_dist ? init_hill_lognor(Y_LN, X, prior):
    init_hill_nor(Y_N, X, prior); 
    break; 
  case cont_model::exp_3:
  case cont_model::exp_5:
    retVal = distribution::log_normal==data_dist ? init_exp_lognor(Y_LN, X, prior):
                                                   init_exp_nor(Y_N, X, prior); 
    break; 
  case cont_model::power: 

    retVal = init_pow_nor( Y_N,  X,  prior);
    break; 
  case cont_model::polynomial:
    /*initialize at Least Squares inv(X.t()*X)*X.t()*Y)
     * 
     */

    deg = distribution::normal_ncv == data_dist? prior.rows() - 3: prior.rows() - 2;   
    retVal = distribution::log_normal==data_dist ? init_poly(Y_LN, X, prior,deg):
                                                   init_poly(Y_N, X, prior,deg); 
    break; 
  default: 
    // this is a different model that shouldn't have MLE fits
    // so we don't do anything
    break; 
  }

  return retVal.col(1);  
}


double compute_lognormal_dof(Eigen::MatrixXd Y,Eigen::MatrixXd X, Eigen::MatrixXd estimate, 
                             bool is_increasing, bool suff_stat, Eigen::MatrixXd prior, 
                             cont_model CM){
  double DOF = 0; 
  Eigen::MatrixXd Xd; 
  Eigen::MatrixXd cv_t; 
  Eigen::MatrixXd pr; 
  Eigen::MatrixXd temp(X.rows(),3);
  Eigen::MatrixXd subBlock(3,3); 
  Eigen::MatrixXd temp_estimate(estimate.rows() + 1,1); 
  
  switch(CM){
  case cont_model::hill:
    Xd = X_gradient_cont<lognormalHILL_BMD_NC>(estimate,Y,X,suff_stat);
    Xd = Xd.block(0,0,Xd.rows(),4); 
    cv_t = X_cov_cont<lognormalHILL_BMD_NC>(estimate,Y,X,suff_stat); 
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    pr = pr.block(0,0,4,4); 
    
    if( fabs(pr.diagonal().array().sum()) == 0){
      DOF = 4.0; 
    }else{
      pr   = Xd.transpose()*cv_t*Xd + pr; 
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    
    break; 
  case cont_model::exp_3:
    
    temp_estimate << estimate(0,0) , estimate(1,0) , 1.0 , estimate.block(2,0,estimate.rows()-2,1); 
    if (is_increasing){
      Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_UP);
      cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_UP);
    }else{
      Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_DOWN);
      cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_DOWN);
    }
    
    temp << Xd.col(0) , Xd.col(1), Xd.col(3);  
    Xd = temp; 
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    subBlock << pr(0,0), pr(0,1), pr(0,3),
                pr(1,0), pr(1,1), pr(1,3),
                pr(3,0), pr(3,1), pr(3,3);
    
    if( fabs(subBlock.diagonal().array().sum()) ==0){
      DOF = 3; 
    } else{
      pr   = Xd.transpose()*cv_t*Xd + subBlock; 
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    break;
  case cont_model::exp_5: 
    if (is_increasing){
      Xd = X_gradient_cont<lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,NORMAL_EXP5_UP);
      cv_t = X_cov_cont< lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat, NORMAL_EXP5_UP);
    }else{
      Xd = X_gradient_cont<lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,NORMAL_EXP5_DOWN);
      cv_t = X_cov_cont< lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat, NORMAL_EXP5_DOWN);
    }
  
    Xd = Xd.block(0,0,Xd.rows(),4); 
   
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    pr = pr.block(0,0,4,4); 
    if( fabs(pr.diagonal().array().sum()) == 0){
      DOF = 4.0; 
    }else{
      pr   = Xd.transpose()*cv_t*Xd + pr; 
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    break; 
  }

  return DOF; 
}

double compute_normal_dof(Eigen::MatrixXd Y,Eigen::MatrixXd X, Eigen::MatrixXd estimate, 
                          bool is_increasing, bool suff_stat,  bool CV, Eigen::MatrixXd prior,
                          cont_model CM,int degree){
  double DOF = 0; 
  Eigen::MatrixXd Xd; 
  Eigen::MatrixXd cv_t; 
  Eigen::MatrixXd pr; 
  Eigen::MatrixXd temp(X.rows(),3);
  Eigen::MatrixXd subBlock(3,3); 
 
  int offset = CV? 1:2; 
  Eigen::MatrixXd temp_estimate(estimate.rows() + 1,1); 
  Eigen::MatrixXd temp_block(1,1); 

  switch(CM){
  case cont_model::polynomial:
    
    Xd = X_gradient_cont_norm<normalPOLYNOMIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,degree);
    temp_block = Xd.block(0,0,Xd.rows(),estimate.rows() - offset); 
    Xd = temp_block; 
    cv_t = X_cov_cont_norm<normalPOLYNOMIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,degree); 
    
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    
    temp_block = pr.block(0,0,estimate.rows() - offset,estimate.rows() - offset); 
    pr = temp_block; 
    
    if( fabs(pr.diagonal().array().sum()) ==0){
      DOF = pr.diagonal().size(); 
    } else{
      temp_block   = Xd.transpose()*cv_t*Xd + pr; 
      pr = temp_block; 
      temp_block = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      Xd = temp_block; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    break; 
  case cont_model::hill:
    Xd = X_gradient_cont_norm<normalHILL_BMD_NC>(estimate,Y,X,suff_stat,CV);
    temp_block  = Xd.block(0,0,Xd.rows(),4); 
    Xd = temp_block; 
    cv_t = X_cov_cont_norm<normalHILL_BMD_NC>(estimate,Y,X,suff_stat,CV); 
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    temp_block  =    pr.block(0,0,4,4); 
    pr = temp_block; 
    
    if( fabs(pr.diagonal().array().sum()) ==0){
      DOF = 4.0; 
    } else{
      temp_block   = Xd.transpose()*cv_t*Xd + pr; 
      pr = temp_block; 
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    
    break; 
  case cont_model::exp_3:
    
    temp_estimate << estimate(0,0) , estimate(1,0) , 1.0 , estimate.block(2,0,estimate.rows()-2,1); 
    if (is_increasing){
      Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_UP);
      cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_UP);
    }else{
      Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_DOWN);
      cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_DOWN);
    }

    temp << Xd.col(0) , Xd.col(1), Xd.col(3);  
    Xd = temp; 
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    subBlock << pr(0,0), pr(0,1), pr(0,3),
                pr(1,0), pr(1,1), pr(1,3),
                pr(3,0), pr(3,1), pr(3,3);
    
    if( fabs(subBlock.diagonal().array().sum()) ==0){
      DOF = 3; 
    } else{
      pr   = Xd.transpose()*cv_t*Xd + subBlock; 
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    break; 
  case cont_model::exp_5: 
    if (is_increasing){
      Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_UP);
      cv_t = X_cov_cont_norm< normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_UP);
    }else{
      Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_DOWN);
      cv_t = X_cov_cont_norm< normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_DOWN);
    }
    
    temp_block = Xd.block(0,0,Xd.rows(),3); 
    Xd = temp_block; 
    
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    temp_block = pr.block(0,0,3,3); 
    pr = temp_block; 
    if( fabs(pr.diagonal().array().sum()) ==0){
      DOF =4.0; 
    } else{
      temp_block   = Xd.transpose()*cv_t*Xd + pr; 
      pr = temp_block; 
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    
    break; 
  case cont_model::power: 
    Xd = X_gradient_cont_norm<normalPOWER_BMD_NC>(estimate,Y,X,CV,suff_stat);
    cv_t = X_cov_cont_norm<normalPOWER_BMD_NC>(estimate,Y,X,CV,suff_stat);
    
    temp_block = Xd.block(0,0,Xd.rows(),3); 
    Xd = temp_block; 
    pr   =  X_logPrior<IDPrior>(estimate,prior); 
    temp_block = pr.block(0,0,3,3); 
    pr = temp_block; 
    
    if( fabs(pr.diagonal().array().sum()) ==0){
      DOF =3.0; 
    }else{
      temp_block   = Xd.transpose()*cv_t*Xd + pr; 
      pr = temp_block;  
      Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
      DOF =  Xd.diagonal().array().sum(); 
    }
    
    break;
  }   
  
  return DOF + offset; 
  
}
  


bool convertSStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                  Eigen::MatrixXd *SSTAT, Eigen::MatrixXd *SSTAT_LN,
                  Eigen::MatrixXd *UX){
  bool convert = true; 
  
  if (Y.cols() == 1 ){
    // check to see if it can be converted into sufficient statistics4
    int temp = 0; 
    // go through each row to see if there are duplicates
    for (int i = 0; i < X.rows(); i++){
      for (int j = 0 ; j < X.rows(); j++){
        if (X(i,0) == X(j,0)){
          temp++; 
        }
      }
      if( temp == 1){
        convert = false; 
      }
      temp = 0; 
    }
    
    if (convert){
      // we can convert the data
       *SSTAT    = createSuffStat( Y, X, false);
       *SSTAT_LN = createSuffStat(Y,X,true); 
        std::vector<double> uniqueX = unique_list(X );
        Eigen::MatrixXd temp_X(uniqueX.size(),1);
        for (int i = 0; i < uniqueX.size(); i++){
          temp_X(i,0) = uniqueX[i]; 
        }
        *UX = temp_X; 
        return true; 
 
     }
  
  }else{
    *SSTAT    = createSuffStat( Y, X, false);
    *SSTAT_LN = createSuffStat(Y,X,true); 
  }
  
  
    
  return false; 
}
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
  
  matrix.conservativeResize(numRows,numCols);
}


void removeCol(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}


bmd_analysis laplace_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                               Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                               bool is_increasing, 
                               double bmrf,   double bk_prob, 
                               double alpha, double step_size,
                               Eigen::MatrixXd init, bool isFast) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
    
  }
  
  
  IDcontinuousPrior model_prior(prior);
  
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp5U(Y, X, suff_stat,  NORMAL_EXP5_UP);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp3U(Y, X, suff_stat,  NORMAL_EXP3_UP);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp5D(Y, X, suff_stat,  NORMAL_EXP5_DOWN);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp3D(Y, X, suff_stat,  NORMAL_EXP3_DOWN);
  lognormalHILL_BMD_NC likelihood_lnhill(Y, X, suff_stat,  0);
  bmd_analysis a;
  Eigen::MatrixXd Xd; 
  Eigen::MatrixXd cv_t; 
  Eigen::MatrixXd pr; 
  Eigen::MatrixXd mean_m; 
  switch (CM)
  {
  
  case cont_model::hill:
#ifdef R_COMPILATION 
    //cout << "Running Hill Model Log-Normality Assumption." << endl;
#endif
    if (isFast){
      
      a =  bmd_fast_BMD_cont <lognormalHILL_BMD_NC, IDcontinuousPrior>
                                (likelihood_lnhill,  model_prior, fixedB, fixedV,
                                 riskType, bmrf, is_increasing,
                                 step_size,init);
    }else{
      a = bmd_analysis_CNC<lognormalHILL_BMD_NC, IDcontinuousPrior>
                              (likelihood_lnhill,  model_prior, fixedB, fixedV,
                                riskType, bmrf, bk_prob,
                                is_increasing, alpha, step_size,init);
    }
    break; 
  case cont_model::exp_3:
#ifdef R_COMPILATION 
    //cout << "Running Exponential 3 Model Log-Normality Assumption." << endl;
#endif 
    if (is_increasing){
      if (isFast){
        a =  bmd_fast_BMD_cont <lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                (likelihood_lnexp3U,  model_prior, fixedB, fixedV,
                                 riskType, bmrf, is_increasing,
                                 step_size,init);
      }else{
      a =  bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                    (likelihood_lnexp3U,  model_prior, fixedB, fixedV,
                                     riskType, bmrf,bk_prob,
                                     is_increasing, alpha, step_size,init);
      }
      
    }else{
      if (isFast){
        a =  bmd_fast_BMD_cont <lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                          (likelihood_lnexp3D,  model_prior, fixedB, fixedV,
                                          riskType, bmrf, is_increasing,
                                          step_size,init);
      }else{
        a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                        (likelihood_lnexp3D,  model_prior, fixedB, fixedV,
                                         riskType, bmrf,bk_prob,
                                         is_increasing, alpha, step_size,init);
      }
    }
    removeRow(a.MAP_ESTIMATE ,2);
    removeRow(a.COV,2);
    removeCol(a.COV,2);
    break; 
  case cont_model::exp_5:
  default: 
#ifdef R_COMPILATION 
    //cout << "Running Exponential 5 Model Log-Normality Assumption." << endl;
#endif
    if (is_increasing){
      if (isFast){
        a =  bmd_fast_BMD_cont <lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                  (likelihood_lnexp5U,  model_prior, fixedB, fixedV,
                                   riskType, bmrf, is_increasing,
                                   step_size,init);
      }else{
        a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                   (likelihood_lnexp5U,  model_prior, fixedB, fixedV,
                                    riskType, bmrf,bk_prob,
                                    is_increasing, alpha, step_size,init);
      }
    }else{
      if (isFast){
        
        a =  bmd_fast_BMD_cont <lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                     (likelihood_lnexp5D,  model_prior, fixedB, fixedV,
                                     riskType, bmrf, is_increasing,
                                     step_size,init);
      }else{
        a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                    (likelihood_lnexp5D,  model_prior, fixedB, fixedV,
                                    riskType, bmrf,bk_prob,
                                    is_increasing, alpha, step_size,init);
      }
    }
  break; 
  
  }
  
  return a; 
}


bmd_analysis laplace_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                            Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                            bool is_increasing, bool bConstVar,
                            double bmrf,   double bk_prob, 
                            double alpha, double step_size, Eigen::MatrixXd init,
                            int degree, bool isFast) {
 
  bool suff_stat = Y.cols() == 1? false:true; 
  
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  
  IDcontinuousPrior model_prior(prior);
  normalPOLYNOMIAL_BMD_NC  likelihood_npoly(Y, X, suff_stat, bConstVar, degree);
  normalHILL_BMD_NC  likelihood_nhill(Y, X, suff_stat, bConstVar, 0);
  normalPOWER_BMD_NC likelihood_power(Y, X, suff_stat, bConstVar, 0);
  normalFUNL_BMD_NC  likelihood_funl(Y, X, suff_stat, bConstVar, 0);
  normalEXPONENTIAL_BMD_NC likelihood_nexp5U(Y, X, suff_stat, bConstVar, NORMAL_EXP5_UP);
  normalEXPONENTIAL_BMD_NC likelihood_nexp3U(Y, X, suff_stat, bConstVar, NORMAL_EXP3_UP);
  normalEXPONENTIAL_BMD_NC likelihood_nexp5D(Y, X, suff_stat, bConstVar, NORMAL_EXP5_DOWN);
  normalEXPONENTIAL_BMD_NC likelihood_nexp3D(Y, X, suff_stat, bConstVar, NORMAL_EXP3_DOWN);

  Eigen::MatrixXd Xd; 
  Eigen::MatrixXd mean_m; 
  Eigen::MatrixXd cv_t; 
  Eigen::MatrixXd pr; 
  bmd_analysis a;
  switch (CM)
  {
    
  case cont_model::funl:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running FUNL Model Normality Assumption using Laplace." << endl;
    }else{
      //cout << "Running FUNL Model Normality-NCV Assumption using Laplace." << endl;
    } 
#endif
    if (isFast){  
      a =  bmd_fast_BMD_cont <normalFUNL_BMD_NC, IDcontinuousPrior>
                            (likelihood_funl,  model_prior, fixedB, fixedV,
                             riskType, bmrf, bk_prob,
                             is_increasing,init);
    }else{
      a = bmd_analysis_CNC<normalFUNL_BMD_NC, IDcontinuousPrior>
                            (likelihood_funl,  model_prior, fixedB, fixedV,
                             riskType, bmrf, bk_prob,
                             is_increasing, alpha, step_size,init);
    }  
  
    break; 
  case cont_model::power:
#ifdef R_COMPILATION      
       if (bConstVar){
            //cout << "Running Power Model Normality Assumption using Laplace." << endl;
       }else{
            //cout << "Running Power Model Normality-NCV Assumption using Laplace." << endl;
       }
#endif
    if (isFast){   
     a =  bmd_fast_BMD_cont <normalPOWER_BMD_NC, IDcontinuousPrior>
                           (likelihood_power,  model_prior, fixedB, fixedV,
                            riskType, bmrf, bk_prob,
                            is_increasing,init);
    }else{
      a = bmd_analysis_CNC<normalPOWER_BMD_NC, IDcontinuousPrior>
                          (likelihood_power,  model_prior, fixedB, fixedV,
                           riskType, bmrf, bk_prob,
                           is_increasing, alpha, step_size,init);
    }
      
    break; 
  case cont_model::hill:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Hill Model Normality Assumption using Laplace." << endl;
    }else{
      //cout << "Running Hill Model Normality-NCV Assumption using Laplace." << endl;
    }
#endif
    if (isFast){
      
        a = bmd_fast_BMD_cont <normalHILL_BMD_NC, IDcontinuousPrior>  (likelihood_nhill,  model_prior, 
                                                             fixedB, fixedV,
                                                             riskType, bmrf, bk_prob,
                                                             is_increasing,init);
    }else{
    a = bmd_analysis_CNC<normalHILL_BMD_NC, IDcontinuousPrior>
                    (likelihood_nhill,  model_prior, fixedB, fixedV,
                     riskType, bmrf, bk_prob,
                     is_increasing, alpha, step_size,init);
    }
    break; 
  case cont_model::exp_3:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Exponential 3 Model Normality Assumption using Laplace." << endl;
    }else{
      //cout << "Running Exponential 3 Model Normality-NCV Assumption using Laplace." << endl;
    }
#endif
    if (is_increasing){
      if (isFast){
        a = bmd_fast_BMD_cont <normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                                                       (likelihood_nexp3U,   model_prior, 
                                                                       fixedB, fixedV,
                                                                       riskType, bmrf, bk_prob,
                                                                       is_increasing,init);
        
        
      }else{
          a =  bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                              (likelihood_nexp3U,  model_prior, fixedB, fixedV,
                               riskType, bmrf,bk_prob,
                               is_increasing, alpha, step_size,init);
      }
    }else{
            if (isFast){
              
              a = bmd_fast_BMD_cont <normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                                          (likelihood_nexp3D,   model_prior, 
                                                           fixedB, fixedV,
                                                           riskType, bmrf, bk_prob,
                                                           is_increasing,init);
                                                       
            }else{
              a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                    (likelihood_nexp3D,  model_prior, fixedB, fixedV,
                     riskType, bmrf,bk_prob,
                     is_increasing, alpha, step_size,init);
            }
    }
    removeRow(a.MAP_ESTIMATE ,2);
    removeRow(a.COV,2);
    removeCol(a.COV,2);
    break; 
  case cont_model::exp_5:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Exponential 5 Model Normality Assumption using Laplace." << endl;
    }else{
      //cout << "Running Exponential 5 Model Normality-NCV Assumption using Laplace." << endl;
    }
#endif
    
    if (is_increasing){
      if (isFast){
            a = bmd_fast_BMD_cont <normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                                              (likelihood_nexp5U,    model_prior, 
                                                               fixedB, fixedV,
                                                               riskType, bmrf, bk_prob,
                                                               is_increasing,init);
            
      }else{
            a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                           (likelihood_nexp5U,  model_prior, fixedB, fixedV,
                                           riskType, bmrf,bk_prob,
                                           is_increasing, alpha, step_size,init);
      }
    }else{
      if (isFast){
          a = bmd_fast_BMD_cont <normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                             (likelihood_nexp5D,   model_prior, 
                                             fixedB, fixedV,
                                             riskType, bmrf, bk_prob,
                                             is_increasing,init);
        
      }else{
            a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                          (likelihood_nexp5D,  model_prior, fixedB, fixedV,
                                          riskType, bmrf,bk_prob,
                                          is_increasing, alpha, step_size,init);
      }
    }
    break; 
  case cont_model::polynomial:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Polynomial Model Normality Assumption using Laplace." << endl;
    }else{
      //cout << "Running Polynomial Model Normality-NCV Assumption using Laplace." << endl;
    }
#endif
    if (isFast){
      a = bmd_fast_BMD_cont <normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>
                                           (likelihood_npoly,   model_prior, 
                                           fixedB, fixedV,
                                           riskType, bmrf, bk_prob,
                                           is_increasing,init);
      
    }else{
    a = bmd_analysis_CNC<normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>
                                          (likelihood_npoly,  model_prior, fixedB, fixedV,
                                           riskType, bmrf, bk_prob,
                                           is_increasing, alpha, step_size,init);
    }
    break;
    
  }
  
  return a; 
}

void transfer_continuous_model(bmd_analysis a, continuous_model_result *model){
	if (model){
	  model->nparms = a.COV.rows(); 
		model->max = a.MAP; 
		model->bmd = a.MAP_BMD; 
		for (int i = 0; i< model->dist_numE; i ++){
				double temp = double(i)/double(model->dist_numE); 
				model->bmd_dist[i] = a.BMD_CDF.inv(temp);     // BMD @ probability
				model->bmd_dist[model->dist_numE + i] = temp; // probability 
		}
		for (int i = 0; i < model->nparms; i++){
				model->parms[i] = a.MAP_ESTIMATE(i,0); 
				for (int j = 0; j < model->nparms; j++){
					model->cov[i + j*model->nparms] = a.COV(i,j); 
				}
		}
	}
}

void inverse_transform_dose(continuous_model_result *model){
  if (model){
    model->bmd = sinh(model->bmd); 
    for (int i = 0; i< model->dist_numE; i ++){
      double temp = double(i)/double(model->dist_numE); 
      model->bmd_dist[i] = sinh(model->bmd_dist[i]);     // BMD @ probability
    }
  }
}


void bmd_range_find(continuousMA_result *res, 
					double *range){
 // assume the minimum BMD for the MA is always 0
 range[0] = 0.0; 
 double current_max = 0.0; 
 for (int j = 10; j > 1; j--){
	 for (int i = 0; i < res->nmodels;i++){
		int temp_idx = res->models[i]->dist_numE -j; 
		
		// make sure we are not dealing with an infinite value
		// or not a number
		if (isfinite(res->models[i]->bmd_dist[temp_idx]) && 
			  !isnan(res->models[i]->bmd_dist[temp_idx])){
			if ( res->models[i]->bmd_dist[temp_idx] > current_max){
				current_max = res->models[i]->bmd_dist[temp_idx]; 
			}
		}
		 
	 }
 }
 // if we don't find a max then the upper limit is NAN
 range[1] = current_max == 0.0 ? std::numeric_limits<double>::quiet_NaN():current_max; 
  
}

void estimate_ma_laplace(continuousMA_analysis *MA,
                         continuous_analysis *CA ,
                         continuousMA_result *res){
  // standardize the data
  
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  bool  tempsa = CA->suff_stat;
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 

  // copy the origional data
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  
  if(!CA->suff_stat){
    //convert to sufficient statistics for speed if we can
    CA->suff_stat = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
    if (CA->suff_stat)// it can be converted
    {
      X = UX; 
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      orig_X = UX;  
      orig_Y = SSTAT; 
      orig_Y_LN = SSTAT_LN;
      
    }else{
      Y = (1/divisor)*Y; // scale the data with the divisor term. 
      X = X/max_dose;
      Y_N = Y; 
      Y_LN = Y; 
    }
  }else{
    orig_Y = cleanSuffStat(Y,X,false,false); 
    orig_Y_LN = cleanSuffStat(Y,X,true,false);
    SSTAT = cleanSuffStat(Y,X,false); 
    SSTAT_LN = cleanSuffStat(Y,X,true);
    
    std::vector<double> tux = unique_list(X); 
    UX = Eigen::MatrixXd(tux.size(),1); 
    for (unsigned int i = 0; i < tux.size(); i++){
      UX(i,0) = tux[i]; 
    }
    Y_N = SSTAT; 
    X = UX; 
    Y_LN = SSTAT_LN; 
  }
  
  
  if (CA->suff_stat){
    X = UX; 

    Eigen::MatrixXd temp; 
    temp = Y_N.col(2);
    Y_N.col(2) = Y_N.col(1);
    Y_N.col(1) = temp; 
    temp = Y_LN.col(2);
    Y_LN.col(2) = Y_LN.col(1);
    Y_LN.col(1) = temp; 
    temp = orig_Y.col(2);
    orig_Y.col(2) = orig_Y.col(1);
    orig_Y.col(1) = temp; 
    temp = orig_Y_LN.col(2);
    orig_Y_LN.col(2) = orig_Y_LN.col(1);
    orig_Y_LN.col(1) = temp; 
    X = X/max_dose;
  } 
  
  bmd_analysis b[MA->nmodels];

#pragma omp parallel
{
  #pragma omp for  
  for (int i = 0; i < MA->nmodels; i++ ){

      std::vector<bool> fixedB; 
      std::vector<double> fixedV;
      // on each iteration make sure there parameters are emptied
      fixedB.clear();
      fixedV.clear(); 
      Eigen::MatrixXd tprior(MA->nparms[i],MA->prior_cols[i]);
      for (int m = 0; m < MA->nparms[i]; m++){
        fixedB.push_back(false);
        fixedV.push_back(0.0); 
        for (int n = 0; n < MA->prior_cols[i]; n++){
          tprior(m,n) = MA->priors[i][m + n*MA->nparms[i]]; 
        }
      }

      Eigen::MatrixXd temp_init =   initialize_model( Y_N, Y_LN, X, 
                                                      tprior,(distribution)MA->disttype[i],(cont_model)MA->models[i]) ;
      
      temp_init = temp_init.array(); 
      
      
      Eigen::MatrixXd init_opt; 
      switch((cont_model)MA->models[i]){
      case cont_model::funl:
         init_opt = bmd_continuous_optimization<normalFUNL_BMD_NC,IDPrior>(Y_N, X, tprior,  fixedB, fixedV, 
                                                                          MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing);
         RescaleContinuousModel<IDPrior>((cont_model)MA->models[i], &tprior, &init_opt, 
                                         1.0, divisor, CA->isIncreasing, MA->disttype[i] == distribution::log_normal,
                                         MA->disttype[i] != distribution::normal_ncv); 
        break; 
      case cont_model::hill:
          init_opt = MA->disttype[i] == distribution::log_normal ?
          bmd_continuous_optimization<lognormalHILL_BMD_NC,IDPrior> (Y_LN, X, tprior, fixedB, fixedV,
                                                                     MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init):
          bmd_continuous_optimization<normalHILL_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                     MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init);
        
          RescaleContinuousModel<IDPrior>((cont_model)MA->models[i], &tprior, &init_opt, 
                                          1.0, divisor, CA->isIncreasing, MA->disttype[i] == distribution::log_normal,
                                          MA->disttype[i] != distribution::normal_ncv); 
          
        break; 
      
      case cont_model::exp_5:
        
 
        init_opt = MA->disttype[i] == distribution::log_normal ?
        bmd_continuous_optimization<lognormalEXPONENTIAL_BMD_NC,IDPrior> (Y_LN, X, tprior, fixedB, fixedV,
                                                                          MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init):
        bmd_continuous_optimization<normalEXPONENTIAL_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                          MA->disttype[i]!= distribution::normal_ncv, CA->isIncreasing,temp_init);

        RescaleContinuousModel<IDPrior>((cont_model)MA->models[i], &tprior, &init_opt, 
                                        1.0, divisor, CA->isIncreasing,MA->disttype[i] == distribution::log_normal,
                                        MA->disttype[i] != distribution::normal_ncv); 
      break; 
      case cont_model::exp_3:
         
          init_opt = MA->disttype[i] == distribution::log_normal ?
        
          bmd_continuous_optimization<lognormalEXPONENTIAL_BMD_NC,IDPrior> (Y_LN, X, tprior, fixedB, fixedV,
                                                                            MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init):
          bmd_continuous_optimization<normalEXPONENTIAL_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                            MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init);
          
          

          RescaleContinuousModel<IDPrior>((cont_model)MA->models[i], &tprior, &init_opt, 
                                          1.0, divisor, CA->isIncreasing,MA->disttype[i] == distribution::log_normal,
                                          MA->disttype[i] != distribution::normal_ncv); 
          
        break; 
      case cont_model::power: 
          init_opt = MA->disttype[i] == distribution::log_normal ?
          bmd_continuous_optimization<lognormalPOWER_BMD_NC,IDPrior> (Y_LN, X, tprior, fixedB, fixedV,
                                                                      MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init):
          bmd_continuous_optimization<normalPOWER_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                      MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing,temp_init);
          
          
          RescaleContinuousModel<IDPrior>((cont_model)MA->models[i], &tprior, &init_opt, 
                                          1.0, divisor, CA->isIncreasing,MA->disttype[i] == distribution::log_normal,
                                          MA->disttype[i] != distribution::normal_ncv); 
          
          break; 
      case cont_model::polynomial:
        init_opt =   bmd_continuous_optimization<normalPOWER_BMD_NC,IDPrior>(Y_N, X, tprior,  fixedB, fixedV, 
                                                                            MA->disttype[i] != distribution::normal_ncv,
                                                                            CA->degree);
        
        RescaleContinuousModel<IDPrior>((cont_model)MA->models[i], &tprior, &init_opt, 
                                        1.0, divisor, CA->isIncreasing,MA->disttype[i] == distribution::log_normal,
                                        MA->disttype[i] != distribution::normal_ncv); 
        break; 
      default:
          break; 
        
      }
     
      // now you fit it based upon the original data
      if (MA->disttype[i] == distribution::log_normal){
        
        
        if (CA->suff_stat ){
          b[i] = laplace_logNormal(orig_Y_LN, X,
                                   tprior, CA->BMD_type, (cont_model)MA->models[i],
                                   CA->isIncreasing, CA->BMR, 
                                   CA->tail_prob,  
                                   CA->alpha, 0.025,init_opt,false);

        }else{
          b[i] = laplace_logNormal(orig_Y_LN, X,
                                   tprior, CA->BMD_type, (cont_model)MA->models[i],
                                   CA->isIncreasing, CA->BMR, 
                                   CA->tail_prob,  
                                   CA->alpha, 0.025,init_opt,false);
          
        }
      
      }else{
        
        bool isNCV = MA->disttype[i] == distribution::normal_ncv? false:true; 
        if (CA->suff_stat ){
           b[i] = laplace_Normal(orig_Y, X,
                                tprior, CA->BMD_type, (cont_model)MA->models[i],
                                CA->isIncreasing,isNCV, CA->BMR, 
                                CA->tail_prob,  
                                CA->alpha, 0.025,init_opt,false);

        }else{
          b[i] = laplace_Normal(orig_Y, X,
                                tprior, CA->BMD_type, (cont_model)MA->models[i],
                                CA->isIncreasing,isNCV, CA->BMR, 
                                CA->tail_prob,  
                                CA->alpha, 0.025,init_opt,false);

        }
        
      }
      
  }
} 

 
  double post_probs[MA->nmodels]; 
  double temp =0.0; 
  double max_prob = -1.0*std::numeric_limits<double>::infinity(); 
  for (int i = 0; i < MA->nmodels; i++){
    temp  = 	b[i].MAP_ESTIMATE.rows()/2 * log(2 * M_PI) - b[i].MAP + 0.5*log(max(0.0,b[i].COV.determinant()));

    if (std::isfinite(temp)){
      max_prob = temp > max_prob? temp:max_prob; 
      post_probs[i] = temp; 
    }else{
      post_probs[i] = -1*std::numeric_limits<double>::infinity();
    }
  }
  
  double norm_sum = 0.0; 

  for (int i = 0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i] - max_prob + log(MA->modelPriors[i]);
    norm_sum     += exp(post_probs[i]);
    post_probs[i] = exp(post_probs[i]);
  }
  
  
  for (int j = 0; j < MA->nmodels; j++){
    
    b[j].COV =  rescale_cov_matrix(b[j].COV, 
                                   b[j].MAP_ESTIMATE, (cont_model) MA->models[j],
                                                                             max_dose, 1.0, false);
    b[j].MAP_ESTIMATE = rescale_parms(b[j].MAP_ESTIMATE,  (cont_model)MA->models[j],
                                      max_dose,1.0, false);
    b[j].MAP_BMD *= max_dose; 
    b[j].BMD_CDF.set_multiple(max_dose);
  
  }
  
  for (int j = 0; j < MA->nmodels; j++){
    post_probs[j] = post_probs[j]/ norm_sum; 
    
    for (double  i = 0.0; i <= 0.50; i += 0.01 ){
      if (!isfinite(b[j].BMD_CDF.inv(i))){
        post_probs[j] = 0;    // if the cdf is infinite before the median
                              // it is removed 
      }  
    } 
  }
  
  norm_sum = 0.0; 
  for (int i =0; i < MA->nmodels; i++){
    norm_sum += post_probs[i]; 
  }
  

  double range[2]; 

  // define the BMD distribution ranges
  // also get compute the MA BMD list
  // this MUST occur after the poster probabilities ABOVE are computed



  for (int i =0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i]/norm_sum; 
    res->post_probs[i] = post_probs[i];
    transfer_continuous_model(b[i],res->models[i]);
    res->models[i]->model = MA->models[i]; 
    res->models[i]->dist  = MA->disttype[i];
  }

  bmd_range_find(res,range);
  
  double prob = 0; 
  int stop = 0; 
  double cbmd = 0.0; 
  do{
       cbmd = (range[1] - range[0])/2; 
       prob = 0; 
       for (int j = 0; j < MA->nmodels; j++){
            if (post_probs[j] > 0){ 
                 prob += b[j].BMD_CDF.P(cbmd)*post_probs[j]; 
            }
       } 
       
       if (prob > 0.99999){
            range[1] = cbmd; 
       }else{
            range[0] = cbmd; 
       } 
       stop ++; 
       
  } while( (prob > 0.99999)  && (stop < 50)); 
  
  range[1] = 2*cbmd; range[0] = 0; 
  
  double trange[2]; 
  trange[0] = 0.0; trange[1] = range[1];
   stop = 0;
  double mid = 0.0; 
  do{
       mid = (trange[1] + trange[0])/2.0; 
       prob = 0; 
       for (int j = 0; j < MA->nmodels; j++){
            if (post_probs[j] > 0){ 
                 prob += b[j].BMD_CDF.P(mid)*post_probs[j]; 
            }
       } 
       
       if (prob < 0.50){
            trange[0] = mid; 
       }else{
            trange[1] = mid;
       }
     
       stop ++; 
       
  } while( (fabs(prob-0.50) > 0.0001) && (stop < 50)); 
   
  double lower_range = mid;
  
  // now find a good starting point for the integration
  mid = 0.0;
  trange[0] = 0.0; trange[1] = range[1];
  do{
       mid = (trange[1] + trange[0])/2.0; 
       prob = 0; 
       for (int j = 0; j < MA->nmodels; j++){
            if (post_probs[j] > 0){ 
                 prob += b[j].BMD_CDF.P(mid)*post_probs[j]; 
            }
       } 
       
       if (prob < 1e-4){
            trange[0] = mid; 
       }else{
            trange[1] = mid;
       }
       
       stop ++; 
       
  } while( (log(fabs(prob-1e-4)) > log(1e-8)) && (stop < 50)); 
  double lower_end = mid; 
  
  
  int i = 0; 
  for (; i < res->dist_numE/2; i ++){
       cbmd =  double(i+1)/double(res->dist_numE/2)*(lower_range-lower_end) + lower_end; 
       double prob = 0.0; 
       
       for (int j = 0; j < MA->nmodels; j++){
            if (post_probs[j] > 0){ 
                 prob += b[j].BMD_CDF.P(cbmd)*post_probs[j]; 
            }
       }
       
       res->bmd_dist[i] = cbmd; 
       res->bmd_dist[i+res->dist_numE]  = prob;
  }
  double range_bmd = range[1] - cbmd;  
 
  
  mid = cbmd;
  
  for (; i < res->dist_numE; i ++){
    double cbmd = mid + (double((i -res->dist_numE/2) + 1)/double(res->dist_numE/2))*(range_bmd); 
    double prob = 0.0; 
    
    for (int j = 0; j < MA->nmodels; j++){
      if (post_probs[j] > 0){ 
          prob += b[j].BMD_CDF.P(cbmd)*post_probs[j]; 
      }
    }
    
    res->bmd_dist[i] = cbmd; 
    res->bmd_dist[i+res->dist_numE]  = prob;
  }
  CA->suff_stat = tempsa;
 
  return;  
}


/*############################################################################
 * 
 * 
 * 
 * 
##############################################################################*/
mcmcSamples mcmc_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                            Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                            bool is_increasing, 
                            double bmrf,   double bk_prob, 
                            double alpha,  int samples, int burnin,  
                            Eigen::MatrixXd initV) {
   
  bool suff_stat = Y.cols() == 1? false:true; 
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }

  mcmcSamples a;
  int adverseR; 
  switch (CM)
  {
  case cont_model::hill:
#ifdef R_COMPILATION 
    //cout << "Running Hill Model Log-Normality Assumption using MCMC." << endl;
#endif
     a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalHILL_BMD_NC, IDPriorMCMC>
                                    (Y,  X, prior, fixedB, fixedV, is_increasing,
                                     bk_prob,suff_stat,bmrf, riskType, alpha,
                                     samples,0,initV); 
    break; 
  case cont_model::exp_3:
      adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
#ifdef R_COMPILATION 
      //cout << "Running Exponential 3 Model Log-Normality Assumption using MCMC." << endl;
#endif
      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXPONENTIAL_BMD_NC, IDPriorMCMC>
                                              (Y,  X, prior, fixedB, fixedV, is_increasing,
                                              bk_prob,suff_stat,bmrf, riskType, alpha,
                                              samples,adverseR,initV);
      // remove the third entry
      removeRow(a.map_cov, 2);
      removeCol(a.map_cov, 2);
      removeRow(a.map_estimate, 2);
    break; 
  case cont_model::exp_5:
  default: 
      adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
#ifdef R_COMPILATION 
      //cout << "Running Exponential 5 Model Log-Normality Assumption using MCMC." << endl;
#endif
      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXPONENTIAL_BMD_NC, IDPriorMCMC>
                                              (Y,  X, prior, fixedB, fixedV, is_increasing,
                                               bk_prob,suff_stat,bmrf, riskType, alpha,
                                               samples,adverseR,initV);
    break; 
    
  }
  
  return a; 
}

mcmcSamples mcmc_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                        Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                        bool is_increasing, bool bConstVar,
                        double bmrf,   double bk_prob, 
                        double alpha, int samples,
                        int burnin, Eigen::MatrixXd initV,
                        int degree = 2); 
mcmcSamples mcmc_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                         Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                         bool is_increasing, bool bConstVar,
                         double bmrf,   double bk_prob, 
                         double alpha, int samples,
                         int burnin, Eigen::MatrixXd initV,
                         int degree) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  
  mcmcSamples a;
  int adverseR; 
  switch (CM)
  {
  case cont_model::funl:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running FUNL Model Normality Assumption using MCMC." << endl;
    }else{
      //cout << "Running FUNL Model Normality-NCV Assumption using MCMC." << endl;
    }
#endif  
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalFUNL_BMD_NC, IDPriorMCMC>
                                            (Y,  X, prior, fixedB, fixedV, is_increasing,
                                             bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                             samples,0,initV); 
    break;   
    
  case cont_model::hill:
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Hill Model Normality Assumption using MCMC." << endl;
    }else{
      //cout << "Running Hill Model Normality-NCV Assumption using MCMC." << endl;
    }
#endif
    
      a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalHILL_BMD_NC, IDPriorMCMC>
                                                (Y,  X, prior, fixedB, fixedV, is_increasing,
                                                 bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                                 samples,0,initV); 
    break; 
  case cont_model::exp_3:
    adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Exponential 3 Model Normality Assumption using MCMC." << endl;
    }else{
      //cout << "Running Exponential 3 Model Normality-NCV Assumption using MCMC." << endl;
    }
#endif
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDPriorMCMC>
                                                (Y,  X, prior, fixedB, fixedV, is_increasing,
                                                 bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                                 samples,adverseR,initV);
    // remove the third entry
    removeRow(a.map_cov, 2);
    removeCol(a.map_cov, 2);
    removeRow(a.map_estimate, 2);
    
    break; 
  case cont_model::exp_5:
    adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Exponential 5 Model Normality Assumption using MCMC." << endl;
    }else{
      //cout << "Running Exponential 5 Model Normality-NCV Assumption using MCMC." << endl;
    }
#endif
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDPriorMCMC>
                                                (Y,  X, prior, fixedB, fixedV, is_increasing,
                                                 bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                                 samples,0,initV);
    break; 
  case cont_model::power:
 
#ifdef R_COMPILATION 
    if (bConstVar){
      //cout << "Running Power Model Normality Assumption using MCMC." << endl;
    }else{
      //cout << "Running Power Model Normality-NCV Assumption using MCMC." << endl;
    }
#endif
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPOWER_BMD_NC, IDPriorMCMC>
                                              (Y,  X, prior, fixedB, fixedV, is_increasing,
                                              bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                              samples,adverseR,initV);
    
  
  break; 
  case cont_model::polynomial:
#ifdef R_COMPILATION 
  if (bConstVar){
    //cout << "Running Polynomial Model Normality Assumption using MCMC." << endl;
  }else{
    //cout << "Running Polynomial Model Normality-NCV Assumption using MCMC." << endl;
  }
#endif

    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPOLYNOMIAL_BMD_NC, IDPriorMCMC>
          (Y,  X, prior, fixedB, fixedV, is_increasing,
           bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
           samples,degree,initV);
    
  default: 
    break; 
  }
  //convert a stuff
  //
  return a; 
}

bmd_analysis create_bmd_analysis_from_mcmc(unsigned int burnin, mcmcSamples s,double max_d){
  bmd_analysis rV;
  rV.MAP          = s.map; 
  rV.MAP_ESTIMATE = s.map_estimate; 
  rV.COV          = s.map_cov; 
  rV.MAP_BMD      = 0; 
  int total = 0; 
  int bad   = 0; 
 
  std::vector<double> v; 
  for (int i = burnin; i < s.BMD.cols(); i++){
    total ++;
    if ( isfinite(s.BMD(0,i)) && s.BMD(0,i) < 10*max_d){ // always look at 5x max dose tested
	        v.push_back(s.BMD(0,i));   // get rid of the burn in samples
    }else{
      bad++; 
    }
  }
  
  double sum = std::accumulate(v.begin(), v.end(), 0.0); 
//= sum/v.size(); //average of the non-infinite bmds
  std::vector<double>  prob;
  std::vector<double> bmd_q; 
  double inf_prob =  double(bad)/double(total); // bad observations 
  if (v.size() > 0){
    std::sort(v.begin(), v.end());
   for (double k = 0.004; k <= 0.9999; k += 0.005){
    	    prob.push_back(k*(1.0-inf_prob)); 
          int idx = int(k*double(v.size()));
          idx = idx == 0? 0: idx-1; 
    	    bmd_q.push_back(v[idx]);
   }
 
    // fix numerical quantile issues.
    for (int i = 1; i < bmd_q.size(); i++){
      if (bmd_q[i] <= bmd_q[i-1]){
         for (int kk = i; kk <  bmd_q.size(); kk++){
            bmd_q[kk] = bmd_q[kk-1] + 1e-6; 
         }
      } 
 
    }
  
    if (prob.size() > 10 && *min_element(bmd_q.begin(), bmd_q.end())  < 1e8 
                         && bmd_q[0] > 0 ){  
        rV.BMD_CDF = bmd_cdf(prob,bmd_q);
    }
    rV.MAP_BMD = rV.BMD_CDF.inv(0.5/(1.0-inf_prob));
  }
   // approximate median; 
  return rV; 
}

void transfer_mcmc_output(mcmcSamples a, bmd_analysis_MCMC *b){

  if (b){
    b->samples = a.samples.cols(); 
    b->nparms  = a.samples.rows(); 

    for(unsigned int i= 0; i < a.BMD.cols();  i++){
      b->BMDS[i] = a.BMD(0,i); 
      for(unsigned int j = 0; j < a.samples.rows(); j++){
        b->parms[i + j*a.BMD.cols()] = a.samples(j,i); 
      }
    }
  }

}

void inverse_transform_dose(bmd_analysis_MCMC *b){
  
  if (b){
    for(unsigned int i= 0; i < b->samples;  i++){
      b->BMDS[i] = sinh(b->BMDS[i]); 
    }
  }
  
}

void estimate_ma_MCMC(continuousMA_analysis *MA,
                      continuous_analysis   *CA,
                      continuousMA_result   *res,
                      ma_MCMCfits           *ma){ 
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  bool tempsa = CA->suff_stat; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 

 
 Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
 Eigen::MatrixXd orig_X = X; 
 
 Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
 Eigen::MatrixXd Y_LN, Y_N;
 
 if(!CA->suff_stat){
   //convert to sufficient statistics for speed if we can
   CA->suff_stat = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
   if (CA->suff_stat)// it can be converted
   {
     X = UX; 
     Y_N = cleanSuffStat(SSTAT,UX,false);  
     Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
     orig_X = UX;  
     orig_Y = SSTAT; 
     orig_Y_LN = SSTAT_LN;
     
   }else{
     Y = (1/divisor)*Y; // scale the data with the divisor term. 
     X = X/max_dose;
     Y_N = Y; 
     Y_LN = Y; 
   }
 }else{
   orig_Y = cleanSuffStat(Y,X,false,false); 
   orig_Y_LN = cleanSuffStat(Y,X,true,false);
   SSTAT = cleanSuffStat(Y,X,false); 
   SSTAT_LN = cleanSuffStat(Y,X,true);
   
   
   std::vector<double> tux = unique_list(X); 
   UX = Eigen::MatrixXd(tux.size(),1); 
   for (unsigned int i = 0; i < tux.size(); i++){
     UX(i,0) = tux[i]; 
   }
   Y_N = SSTAT; 
   X = UX; 
   Y_LN = SSTAT_LN; 
 }
 
 
 
 if (CA->suff_stat){
   X = UX; 
   
   Eigen::MatrixXd temp; 
   temp = Y_N.col(2);
   Y_N.col(2) = Y_N.col(1);
   Y_N.col(1) = temp; 
   temp = Y_LN.col(2);
   Y_LN.col(2) = Y_LN.col(1);
   Y_LN.col(1) = temp; 
   temp = orig_Y.col(2);
   orig_Y.col(2) = orig_Y.col(1);
   orig_Y.col(1) = temp; 
   temp = orig_Y_LN.col(2);
   orig_Y_LN.col(2) = orig_Y_LN.col(1);
   orig_Y_LN.col(1) = temp; 
   X = X/max_dose;
 }
 
  
  mcmcSamples a[MA->nmodels];

  unsigned int samples = CA->samples; 
  unsigned int burnin  = CA->burnin;  
  
#pragma omp parallel
{
#pragma omp for
  for (int i = 0; i < MA->nmodels; i++ ){
    std::vector<bool> fixedB; 
    std::vector<double> fixedV; 
    fixedB.clear(); // on each iteration make sure there parameters are emptied
    fixedV.clear(); 
    Eigen::MatrixXd tprior(MA->nparms[i],MA->prior_cols[i]);
    for (int m = 0; m < MA->nparms[i]; m++){
        fixedB.push_back(false);
        fixedV.push_back(0.0); 
        for (int n = 0; n < MA->prior_cols[i]; n++){
          tprior(m,n) = MA->priors[i][m + n*MA->nparms[i]]; 
        }
    }
  Eigen::MatrixXd temp_init =   initialize_model( Y_N, Y_LN, X, 
                                                    tprior,(distribution)MA->disttype[i],(cont_model)MA->models[i]) ;
    
   temp_init = temp_init.array(); 
    
   Eigen::MatrixXd init_opt; 
   switch((cont_model)MA->models[i]){
     case cont_model::funl:
       init_opt =  bmd_continuous_optimization<normalFUNL_BMD_NC,IDPriorMCMC> (Y_N, X, tprior,  fixedB, fixedV, 
                                                                           MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing);
       
       RescaleContinuousModel<IDPriorMCMC>((cont_model)MA->models[i], &tprior, &init_opt, 
                                       1.0, divisor, CA->isIncreasing, MA->disttype[i] == distribution::log_normal,
                                       MA->disttype[i] != distribution::normal_ncv); 
       
       break; 
       case cont_model::hill:
          init_opt = MA->disttype[i] == distribution::log_normal ?
                     bmd_continuous_optimization<lognormalHILL_BMD_NC,IDPriorMCMC> (Y_LN, X, tprior, fixedB, fixedV,
                                                                                MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing, temp_init):
                     bmd_continuous_optimization<normalHILL_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                                MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing, temp_init);
          //updated prior updated 
           RescaleContinuousModel<IDPriorMCMC>((cont_model)MA->models[i], &tprior, &init_opt, 
                                          1.0, divisor, CA->isIncreasing, MA->disttype[i] == distribution::log_normal,
                                          MA->disttype[i] != distribution::normal_ncv); 
           
         break; 
       case cont_model::exp_3:
       case cont_model::exp_5:
         init_opt = MA->disttype[i] == distribution::log_normal ?
                   bmd_continuous_optimization<lognormalEXPONENTIAL_BMD_NC,IDPriorMCMC> (Y_LN, X, tprior, fixedB, fixedV,
                                                                                     MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing, temp_init):
                   bmd_continuous_optimization<normalEXPONENTIAL_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                                     MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing, temp_init);
         //updated prior updated 
         //updated prior updated 
         RescaleContinuousModel<IDPriorMCMC>((cont_model) MA->models[i], &tprior, &init_opt, 
                                         1.0, divisor, CA->isIncreasing,MA->disttype[i] == distribution::log_normal,
                                         MA->disttype[i] != distribution::normal_ncv); 
                   
         break; 
       case cont_model::power: 
         init_opt = MA->disttype[i] == distribution::log_normal ?
                   bmd_continuous_optimization<lognormalPOWER_BMD_NC,IDPriorMCMC> (Y_LN, X, tprior, fixedB, fixedV,
                                                                               MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing, temp_init):
                   bmd_continuous_optimization<normalPOWER_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                               MA->disttype[i] != distribution::normal_ncv, CA->isIncreasing, temp_init);
                 
         //updated prior updated 
         RescaleContinuousModel<IDPriorMCMC>((cont_model)MA->models[i], &tprior, &init_opt, 
                                                   1.0, divisor, CA->isIncreasing,MA->disttype[i] == distribution::log_normal,
                                                   MA->disttype[i] != distribution::normal_ncv); 
         
         break; 
       case cont_model::polynomial:
       default:
           break; 
         
     }
      
    a[i] = MA->disttype[i] == distribution::log_normal?
                              mcmc_logNormal(orig_Y_LN, X,
                                                tprior, CA->BMD_type, (cont_model)MA->models[i],
                                                CA->isIncreasing, CA->BMR, 
                                                CA->tail_prob,  
                                                CA->alpha, samples, burnin,init_opt):
                              mcmc_Normal(orig_Y, X,
                                              tprior, CA->BMD_type, (cont_model)MA->models[i],
                                              CA->isIncreasing, MA->disttype[i] != distribution::normal_ncv, CA->BMR,  
                                              CA->tail_prob,  
                                              CA->alpha,samples, burnin,init_opt);
         
      
  }
}  

  bmd_analysis b[MA->nmodels]; 
  double temp_m_dose = orig_X.maxCoeff();
  for (int i = 0; i < MA->nmodels; i++){
    b[i] = create_bmd_analysis_from_mcmc(burnin,a[i],temp_m_dose);
  }

  double post_probs[MA->nmodels]; 
  double temp =0.0; 
  double max_prob = -1.0*std::numeric_limits<double>::infinity(); 
  for (int i = 0; i < MA->nmodels; i++){
    temp  = 	b[i].MAP_ESTIMATE.rows()/2 * log(2 * M_PI) - b[i].MAP + 0.5*log(max(0.0,b[i].COV.determinant()));
    if (std::isfinite(temp)){
      max_prob = temp > max_prob? temp:max_prob; 
      post_probs[i] = temp; 
    }else{
      post_probs[i] = -1*std::numeric_limits<double>::infinity();
    } 
    post_probs[i] = temp; 
  }
  
  double norm_sum = 0.0; 
 
  for (int i = 0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i] - max_prob + log(MA->modelPriors[i]); //FIX ME: ADD MODEL PROBS
    norm_sum     += exp(post_probs[i]);
    post_probs[i] = exp(post_probs[i]);
  }

  for (int j = 0; j < MA->nmodels; j++){
    post_probs[j] = post_probs[j]/ norm_sum; 

    for (double  i = 0.0; i <= 0.5; i += 0.01 ){
      if ( !isfinite(b[j].BMD_CDF.inv(i)) || isnan(b[j].BMD_CDF.inv(i))){
        
         post_probs[j] = 0;    // if the cdf has nan/inf before the median
                               // it is removed from the calculation and given a 0 posterior
      }  
    } 
  }
  
  norm_sum = 0.0; 
  for (int i =0; i < MA->nmodels; i++){
    norm_sum += post_probs[i]; 
  }
  
  
  /////////////////////////////////////////
  /////////// Rescale mcmc
  for (int  i = 0; i < MA->nmodels; i++){
      rescale_mcmc(&a[i], (cont_model)MA->models[i],
                   max_dose, MA->disttype[i] , 2);
  }
  
  
  for (int j = 0; j < MA->nmodels; j++){
    b[j].COV =  rescale_cov_matrix(b[j].COV, 
                                   b[j].MAP_ESTIMATE, (cont_model)MA->models[j],
                                                                            max_dose, 1.0, false);
    
    b[j].MAP_ESTIMATE = rescale_parms(b[j].MAP_ESTIMATE,  (cont_model)MA->models[j],
                                      max_dose,1.0, false);
    b[j].MAP_BMD *= max_dose; 
    b[j].BMD_CDF.set_multiple(max_dose);
  }
   ////////////////////////////////////////////
  for (int i =0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i]/norm_sum; 
    res->post_probs[i] = post_probs[i];
    transfer_continuous_model(b[i],res->models[i]);
    transfer_mcmc_output(a[i],ma->analyses[i]); 
    res->models[i]->model = MA->models[i]; 
    res->models[i]->dist  = MA->disttype[i];
   
  }
 
  
  double range[2]; 
  
  // define the BMD distribution ranges
  // also get compute the MA BMD list
  bmd_range_find(res,range);
 
 double prob = 0; 
 int stop = 0; 
 double cbmd = 0.0; 
 do{
      cbmd = (range[1] - range[0])/2; 
      prob = 0; 
      for (int j = 0; j < MA->nmodels; j++){
           if (post_probs[j] > 0){ 
                prob += b[j].BMD_CDF.P(cbmd)*post_probs[j]; 
           }
      } 
      
      if (prob > 0.99999){
           range[1] = cbmd; 
      }else{
           range[0] = cbmd; 
      } 
      stop ++; 
      
 } while( (prob > 0.99999)  && (stop < 16)); 
 
 range[1] = 2*cbmd; range[0] = 0; 
 
 double trange[2]; 
 trange[0] = 0.0; trange[1] = range[1];
 stop = 0;
 double mid = 0.0; 
 do{
      mid = (trange[1] + trange[0])/2.0; 
      prob = 0; 
      for (int j = 0; j < MA->nmodels; j++){
           if (post_probs[j] > 0){ 
                prob += b[j].BMD_CDF.P(mid)*post_probs[j]; 
           }
      } 
      
      if (prob < 0.50){
           trange[0] = mid; 
      }else{
           trange[1] = mid;
      }
      
      stop ++; 
      
 } while( (fabs(prob-0.50) > 0.0001) && (stop < 50)); 
 
 double lower_range = mid;
 // now find a good starting point for the integration
 mid = 0.0;
 trange[0] = 0.0; trange[1] = range[1];
 do{
      mid = (trange[1] + trange[0])/2.0; 
      prob = 0; 
      for (int j = 0; j < MA->nmodels; j++){
           if (post_probs[j] > 0){ 
                prob += b[j].BMD_CDF.P(mid)*post_probs[j]; 
           }
      } 
      
      if (prob < 1e-4){
           trange[0] = mid; 
      }else{
           trange[1] = mid;
      }
      
      stop ++; 
      
 } while( (log(fabs(prob-1e-4)) > log(1e-8)) && (stop < 50)); 
 double lower_end = mid; 
 
 int i = 0; 
 for (; i < res->dist_numE/2; i ++){
      cbmd =  double(i+1)/double(res->dist_numE/2)*(lower_range-lower_end)+lower_end; 
      double prob = 0.0; 
      
      for (int j = 0; j < MA->nmodels; j++){
           if (post_probs[j] > 0){ 
                prob += b[j].BMD_CDF.P(cbmd)*post_probs[j]; 
           }
      }
      
      res->bmd_dist[i] = cbmd; 
      res->bmd_dist[i+res->dist_numE]  = prob;
 }
 double range_bmd = range[1] - cbmd;  
 
 
 mid = cbmd;
 
 for (; i < res->dist_numE; i ++){
      double cbmd = mid + (double((i -res->dist_numE/2) + 1)/double(res->dist_numE/2))*(range_bmd); 
      double prob = 0.0; 
      
      for (int j = 0; j < MA->nmodels; j++){
           if (post_probs[j] > 0){ 
                prob += b[j].BMD_CDF.P(cbmd)*post_probs[j]; 
           }
      }
      
      res->bmd_dist[i] = cbmd; 
      res->bmd_dist[i+res->dist_numE]  = prob;
 }
 
 CA->suff_stat = tempsa;
  return; 
}

/*estimate a single model using laplace/profile likelihood*/
void estimate_sm_laplace(continuous_analysis *CA ,
                         continuous_model_result *res, bool isFast){
  // standardize the data
  
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  bool tempsa = CA->suff_stat;
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 

  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->transform_dose?asinh(CA->doses[i]):CA->doses[i]; 
    if(CA->suff_stat){
      
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  

  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  
 

  if(!CA->suff_stat){
    //convert to sufficient statistics for speed if we can
    CA->suff_stat = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
    if (CA->suff_stat)// it can be converted
    {
      X = UX; 
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      orig_X = UX;  
      orig_Y = SSTAT; 
      orig_Y_LN = SSTAT_LN;
      
    }else{
      Y = (1/divisor)*Y; // scale the data with the divisor term. 
      X = X/max_dose;
      Y_N = Y; 
      Y_LN = Y; 
    }
  }else{
    orig_Y = cleanSuffStat(Y,X,false,false); 
    orig_Y_LN = cleanSuffStat(Y,X,true,false);
    SSTAT = cleanSuffStat(Y,X,false); 
    SSTAT_LN = cleanSuffStat(Y,X,true);
    
    
    std::vector<double> tux = unique_list(X); 
    UX = Eigen::MatrixXd(tux.size(),1); 
    for (unsigned int i = 0; i < tux.size(); i++){
      UX(i,0) = tux[i]; 
    }
    Y_N = SSTAT; 
    X = UX; 
    Y_LN = SSTAT_LN; 
  }
  
  if (CA->suff_stat){
    X = UX; 
    
    Eigen::MatrixXd temp; 
    temp = Y_N.col(2);
    Y_N.col(2) = Y_N.col(1);
    Y_N.col(1) = temp; 
    temp = Y_LN.col(2);
    Y_LN.col(2) = Y_LN.col(1);
    Y_LN.col(1) = temp; 
    temp = orig_Y.col(2);
    orig_Y.col(2) = orig_Y.col(1);
    orig_Y.col(1) = temp; 
    temp = orig_Y_LN.col(2);
    orig_Y_LN.col(2) = orig_Y_LN.col(1);
    orig_Y_LN.col(1) = temp; 
    X = X/max_dose;
  }
 
  
  bmd_analysis b; 
  std::vector<bool> fixedB; 
  std::vector<double> fixedV; 

  Eigen::MatrixXd tprior(CA->parms,CA->prior_cols);
  for (int m = 0; m < CA->parms; m++){
    fixedB.push_back(false);
    fixedV.push_back(0.0); 
    for (int n = 0; n < CA->prior_cols; n++){
      tprior(m,n) = CA->prior[m + n*CA->parms]; 
    }
  }
  
  Eigen::MatrixXd temp_init =   initialize_model( Y_N, Y_LN, X, 
                                                  tprior,(distribution)CA->disttype,CA->model) ;
 
  
  temp_init = temp_init.array(); 
  
  Eigen::MatrixXd init_opt; 

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  
  /* Getting number of milliseconds as an integer. */
  
  switch((cont_model)CA->model){
  case cont_model::funl:
 
    init_opt =  bmd_continuous_optimization<normalFUNL_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                           CA->disttype != distribution::normal_ncv,
                                                                           CA->isIncreasing,temp_init);

    RescaleContinuousModel<IDPrior>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, CA->isIncreasing, CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 

  
    break; 
    
  case cont_model::hill:

    if( CA->disttype == distribution::log_normal ){
      init_opt = bmd_continuous_optimization<lognormalHILL_BMD_NC,IDPrior> (Y_LN, X, tprior, fixedB, fixedV,
                                                              CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init);
    }else{
      init_opt = bmd_continuous_optimization<normalHILL_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                               CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init);

    
    }
    RescaleContinuousModel<IDPrior>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, CA->isIncreasing, CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
    
    
    break; 
  case cont_model::exp_3:
  case cont_model::exp_5:
  
    init_opt = CA->disttype == distribution::log_normal ?
    bmd_continuous_optimization<lognormalEXPONENTIAL_BMD_NC,IDPrior> (Y_LN, X, tprior, fixedB, fixedV,
                                                                      CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init):
    bmd_continuous_optimization<normalEXPONENTIAL_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                      CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init);
    RescaleContinuousModel<IDPrior>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, CA->isIncreasing, CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
    
    
    break; 
  case cont_model::power: 
    init_opt =  bmd_continuous_optimization<normalPOWER_BMD_NC,IDPrior>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                           CA->disttype != distribution::normal_ncv, 
                                                                           CA->isIncreasing,temp_init);
   
    RescaleContinuousModel<IDPrior>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, 
                                    CA->isIncreasing,CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
    
     
    break; 
  case cont_model::polynomial:
    // Polynomials are ONLY normal models. 
    
    init_opt =  bmd_continuous_optimization<normalPOLYNOMIAL_BMD_NC,IDPrior> (Y_N, X, tprior,  fixedB, fixedV, 
                                                                         CA->disttype != distribution::normal_ncv,
                                                                         CA->isIncreasing,
                                                                         CA->parms - 2 - (CA->disttype == distribution::normal_ncv ));

    RescaleContinuousModel<IDPrior>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, 
                                    CA->isIncreasing,CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
  default:
    break; 
    
  }
 
 
  double DOF = 0; 
  // now you fit it based upon the origional data
  
  if (CA->disttype == distribution::log_normal){
     if (CA->suff_stat ){
      
      b = laplace_logNormal(orig_Y_LN, X,
                            tprior, CA->BMD_type, (cont_model)CA->model,
                            CA->isIncreasing, CA->BMR, 
                            CA->tail_prob,  
                            CA->alpha, 0.025,init_opt,isFast);
    }else{
      
      b = laplace_logNormal(orig_Y_LN, X,
                            tprior, CA->BMD_type, (cont_model)CA->model,
                            CA->isIncreasing, CA->BMR, 
                            CA->tail_prob,  
                            CA->alpha, 0.025,init_opt,isFast);
      
    }

     DOF =  compute_lognormal_dof(orig_Y_LN,X, b.MAP_ESTIMATE, 
                                 CA->isIncreasing, CA->suff_stat, tprior, 
                                 CA->model); 
    
  }else{

    bool isNCV = CA->disttype != distribution::normal_ncv; 
     if (CA->suff_stat ){
       b = laplace_Normal(orig_Y, X,
                        tprior, CA->BMD_type, (cont_model) CA->model,
                        CA->isIncreasing,isNCV, CA->BMR, 
                        CA->tail_prob,  
                        CA->alpha, 0.025,init_opt,CA->degree,isFast);
      
    }else{

      b = laplace_Normal(orig_Y, X,
                         tprior, CA->BMD_type, (cont_model)CA->model,
                         CA->isIncreasing,isNCV, CA->BMR, 
                         CA->tail_prob,  
                         CA->alpha, 0.025,init_opt,CA->degree,isFast);
        
    }
    
   DOF =  compute_normal_dof(orig_Y,X, b.MAP_ESTIMATE, 
                             CA->isIncreasing, CA->suff_stat, isNCV,tprior, 
                             CA->model,CA->degree); 
    
    
  }
  

   ///////////////////////////////////////////////////////
  // NOTE: need to rescale the covariance first

  b.COV =  rescale_cov_matrix(b.COV, 
                              b.MAP_ESTIMATE, CA->model,
                              max_dose, 1.0, false,CA->degree);

  b.MAP_ESTIMATE = rescale_parms(b.MAP_ESTIMATE, CA->model,
                                 max_dose,1.0, false, CA->degree);
  b.MAP_BMD *= max_dose; 
  b.BMD_CDF.set_multiple(max_dose); 
  ///////////////////////////////////////////////////////
  std::vector<double> v(orig_X.rows()); 
  for (int i ; i < orig_X.rows(); i++){
    v[i] = orig_X(i,0); 
  } 
  transfer_continuous_model(b,res);
  // transform back to the origional scale
  if (CA->transform_dose){
    inverse_transform_dose(res); 
  }
  
  res->model_df = DOF; 
  res->total_df = std::set<double>( v.begin(), v.end() ).size() - DOF; 

  res->model = CA->model; 
  res->dist  = CA->disttype; 
  CA->suff_stat = tempsa;

  return;  
}


void estimate_sm_mcmc(continuous_analysis *CA,
                      continuous_model_result *res,
                      bmd_analysis_MCMC *mcmc)
{ 
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  bool tempsa = CA->suff_stat; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->transform_dose?asinh(CA->doses[i]):CA->doses[i];  
    if(CA->suff_stat){
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }

  
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  
  if(!CA->suff_stat){
    //convert to sufficient statistics for speed if we can
    CA->suff_stat = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
    if (CA->suff_stat)// it can be converted
    {
      X = UX; 
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      orig_X = UX;  
      orig_Y = SSTAT; 
      orig_Y_LN = SSTAT_LN;
      
    }else{
      Y = (1/divisor)*Y; // scale the data with the divisor term. 
      X = X/max_dose;
      Y_N = Y; 
      Y_LN = Y; 
    }
  }else{
    orig_Y = cleanSuffStat(Y,X,false,false); 
    orig_Y_LN = cleanSuffStat(Y,X,true,false);
    SSTAT = cleanSuffStat(Y,X,false); 
    SSTAT_LN = cleanSuffStat(Y,X,true);
    
    std::vector<double> tux = unique_list(X); 
    UX = Eigen::MatrixXd(tux.size(),1); 
    for (unsigned int i = 0; i < tux.size(); i++){
      UX(i,0) = tux[i]; 
    }
    Y_N = SSTAT; 
    X = UX; 
    Y_LN = SSTAT_LN; 
  }

  if (CA->suff_stat){
    X = UX; 
    //  Y_N = cleanSuffStat(SSTAT,UX,false);  
    //  Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
    Eigen::MatrixXd temp; 
    temp = Y_N.col(2);
    Y_N.col(2) = Y_N.col(1);
    Y_N.col(1) = temp; 
    temp = Y_LN.col(2);
    Y_LN.col(2) = Y_LN.col(1);
    Y_LN.col(1) = temp; 
    temp = orig_Y.col(2);
    orig_Y.col(2) = orig_Y.col(1);
    orig_Y.col(1) = temp; 
    temp = orig_Y_LN.col(2);
    orig_Y_LN.col(2) = orig_Y_LN.col(1);
    orig_Y_LN.col(1) = temp; 
    X = X/max_dose;
  }
  
  mcmcSamples a;
  unsigned int samples = CA->samples; 
  unsigned int burnin  = CA->burnin;  
  
  std::vector<bool> fixedB; 
  std::vector<double> fixedV; 
  
  Eigen::MatrixXd tprior(CA->parms,CA->prior_cols);
  for (int m = 0; m < CA->parms; m++){
       fixedB.push_back(false);
       fixedV.push_back(0.0); 
       for (int n = 0; n < CA->prior_cols; n++){
            tprior(m,n) = CA->prior[m + n*CA->parms]; 
       }
  }
  bool bob = CA->disttype != distribution::normal_ncv; 
 
  Eigen::MatrixXd temp_init =   initialize_model( Y_N, Y_LN, X, 
                                                  tprior,(distribution)CA->disttype,CA->model ) ;
  temp_init = temp_init.array(); 
  
  Eigen::MatrixXd init_opt; 
  
  switch((cont_model)CA->model){
  case cont_model::funl:
    
    init_opt =  bmd_continuous_optimization<normalFUNL_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                           CA->disttype != distribution::normal_ncv,
                                                                           CA->isIncreasing,temp_init);
    
    RescaleContinuousModel<IDPriorMCMC>((cont_model)CA->model, &tprior, &init_opt, 
                                       1.0, divisor, CA->isIncreasing, CA->disttype == distribution::log_normal,
                                        CA->disttype != distribution::normal_ncv); 
    
    
    break; 
  case cont_model::hill:
    init_opt = CA->disttype == distribution::log_normal ?
    bmd_continuous_optimization<lognormalHILL_BMD_NC,IDPriorMCMC> (Y_LN, X, tprior, fixedB, fixedV,
                                                               CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init):
    bmd_continuous_optimization<normalHILL_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                               CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init);
    
    RescaleContinuousModel<IDPriorMCMC>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, CA->isIncreasing, CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
    
    break; 
  case cont_model::exp_3:
  case cont_model::exp_5:
    
    init_opt = CA->disttype == distribution::log_normal ?
    bmd_continuous_optimization<lognormalEXPONENTIAL_BMD_NC,IDPriorMCMC> (Y_LN, X, tprior, fixedB, fixedV,
                                                                      CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init):
    bmd_continuous_optimization<normalEXPONENTIAL_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                      CA->disttype != distribution::normal_ncv, CA->isIncreasing,temp_init);
   
    RescaleContinuousModel<IDPriorMCMC>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, CA->isIncreasing, CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
    
    break; 
  case cont_model::power: 
    init_opt =  bmd_continuous_optimization<normalPOWER_BMD_NC,IDPriorMCMC>    (Y_N, X, tprior,  fixedB, fixedV, 
                                                                            CA->disttype != distribution::normal_ncv, 
                                                                            CA->isIncreasing,temp_init);
    
    RescaleContinuousModel<IDPriorMCMC>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, 
                                    CA->isIncreasing,CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
    
    
    break; 
  case cont_model::polynomial:
    // Polynomials are ONLY normal models. 
    init_opt =  bmd_continuous_optimization<normalPOLYNOMIAL_BMD_NC,IDPriorMCMC> (Y_N, X, tprior,  fixedB, fixedV, 
                                                                              CA->disttype != distribution::normal_ncv,
                                                                              CA->isIncreasing,
                                                                              CA->parms - 2 - (CA->disttype == distribution::normal_ncv ));
    RescaleContinuousModel<IDPriorMCMC>((cont_model)CA->model, &tprior, &init_opt, 
                                    1.0, divisor, 
                                    CA->isIncreasing,CA->disttype == distribution::log_normal,
                                    CA->disttype != distribution::normal_ncv); 
  default:
    break; 
  
  }
 
    
  a = CA->disttype == distribution::log_normal?
       mcmc_logNormal(orig_Y_LN, X,
                     tprior, CA->BMD_type, (cont_model)CA->model,
                     CA->isIncreasing, CA->BMR, 
                     CA->tail_prob,  
                     CA->alpha, samples, burnin,init_opt):
       mcmc_Normal(orig_Y, X,
                    tprior, CA->BMD_type, (cont_model)CA->model,
                    CA->isIncreasing, CA->disttype != distribution::normal_ncv, CA->BMR,  
                    CA->tail_prob,  
                    CA->alpha,samples, burnin,init_opt,CA->degree);
  
  ///////////////////////////////////
  /////////// Rescale mcmc
  rescale_mcmc(&a, (cont_model)CA->model,
                max_dose, CA->disttype,CA->degree);
  ///////////////////////////////////
  

  bmd_analysis b; 
  CA->suff_stat = tempsa;
  b = create_bmd_analysis_from_mcmc(burnin,a,max_dose);
  
  transfer_continuous_model(b,res);
  if (CA->transform_dose){
    inverse_transform_dose(res); 
  }
  
  transfer_mcmc_output(a,mcmc); 
  if (CA->transform_dose){
    inverse_transform_dose(res); 
    inverse_transform_dose(mcmc);
  }
  res->model = CA->model; 
  res->dist  = CA->disttype;
 
  return; 
}


/*
 * 
 */
void estimate_log_normal_aod(continuous_analysis *CA,
                             continuous_deviance *aod){
  
  
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  bool can_be_suff = true; 
  if (Y.cols() == 1){
     can_be_suff = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
  }
  else{
    SSTAT     = cleanSuffStat(Y,X,false,false); 
    SSTAT_LN  = cleanSuffStat(Y,X,true,false);
    UX = X; 
  }
 
  if(!can_be_suff){
     aod->R  =  std::numeric_limits<double>::infinity();
     aod->A1 =  std::numeric_limits<double>::infinity();
     aod->A2 =  std::numeric_limits<double>::infinity();
     aod->A3 =  std::numeric_limits<double>::infinity();
     return;   
  }else{
     Y_LN = SSTAT_LN; 
     Eigen::MatrixXd temp = Y_LN.col(2);
     Y_LN.col(2) = Y_LN.col(1);
     Y_LN.col(1) = temp; 
    log_normal_AOD_fits(Y_LN, UX, 
                        can_be_suff, aod);
    return; 
  }
}


/*
 * 
 */
void estimate_normal_aod(continuous_analysis *CA,
                         continuous_deviance *aod){
  
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  bool can_be_suff = true; 
  if (Y.cols() == 1){ 
    //individual data
     can_be_suff = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
  }else{
    
    SSTAT     = cleanSuffStat(Y,X,false,false); 
    SSTAT_LN  = cleanSuffStat(Y,X,true,false);
    UX = X; 
  }
 
  
  if(!can_be_suff){
    
    
    aod->R  =  std::numeric_limits<double>::infinity();
    aod->A1 =  std::numeric_limits<double>::infinity();
    aod->A2 =  std::numeric_limits<double>::infinity();
    aod->A3 =  std::numeric_limits<double>::infinity();
  
    return;   
  }else{
    Y_N = SSTAT; 
    Eigen::MatrixXd temp = Y_N.col(2);
    Y_N.col(2) = Y_N.col(1);
    Y_N.col(1) = temp; 
    normal_AOD_fits(Y_N, UX, 
                    can_be_suff, aod);
    return; 
  }
}

void estimate_normal_variance(continuous_analysis *CA,
                              double *v_c, double *v_nc, double *v_pow){
  
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  bool can_be_suff = true; 
  if (Y.cols() == 1){ 
    //individual data
    can_be_suff = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
  }else{
    
    SSTAT     = cleanSuffStat(Y,X,false,false); 
    SSTAT_LN  = cleanSuffStat(Y,X,true,false);
    UX = X; 
  }

  
  if(!can_be_suff){
    *v_c   =  std::numeric_limits<double>::infinity();
    *v_nc  =  std::numeric_limits<double>::infinity();
    *v_pow =  std::numeric_limits<double>::infinity();
    return;   
  }else{
     Y_N = SSTAT; 
     Eigen::MatrixXd temp = Y_N.col(2);
     Y_N.col(2) = Y_N.col(1);
     Y_N.col(1) = temp; 
     variance_fits(Y_N, UX, 
                    can_be_suff, 
                    v_c, v_nc, v_pow);
    
    return; 
  }
}

                   
/*
 * 
 * 
 */
void continuous_expectation( const continuous_analysis *CA, const continuous_model_result *MR,
                             continuous_expected_result *expected){

  // copy the origional data
  
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  bool tempsa = CA->suff_stat; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      Y(i,2) = CA->sd[i]; 
      Y(i,1) = CA->n_group[i]; 
    }
  }
  
  Eigen::MatrixXd myX = X; 
  double divisor = get_divisor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
  Eigen::MatrixXd orig_X = X; 
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  bool suff_stat = false; 
  if(!CA->suff_stat){
    //convert to sufficient statistics for speed if we can
    suff_stat  = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
    if (suff_stat)// it can be converted
    {
      X = UX; 
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      orig_X = UX;  
      orig_Y = SSTAT; 
      orig_Y_LN = SSTAT_LN;
      
    }else{
      Y = Y; // scale the data with the divisor term.
      Y_N = Y; 
      Y_LN = Y; 
    }
  }else{
    suff_stat = true; 
    orig_Y = cleanSuffStat(Y,X,false,false); 
    orig_Y_LN = cleanSuffStat(Y,X,true,false);
    SSTAT = cleanSuffStat(Y,X,false); 
    SSTAT_LN = cleanSuffStat(Y,X,true);
    
    std::vector<double> tux = unique_list(X); 
    UX = Eigen::MatrixXd(tux.size(),1); 
    for (unsigned int i = 0; i < tux.size(); i++){
      UX(i,0) = tux[i]; 
    }
    Y_N = SSTAT; 
    X = UX; 
    Y_LN = SSTAT_LN; 
  }
  
  if (suff_stat){
  
    X = UX; 
    //  Y_N = cleanSuffStat(SSTAT,UX,false);  
    //  Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
    Eigen::MatrixXd temp; 
    temp = Y_N.col(2);
    Y_N.col(2) = Y_N.col(1);
    Y_N.col(1) = temp; 
    temp = Y_LN.col(2);
    Y_LN.col(2) = Y_LN.col(1);
    Y_LN.col(1) = temp; 
    temp = orig_Y.col(2);
    orig_Y.col(2) = orig_Y.col(1);
    orig_Y.col(1) = temp; 
    temp = orig_Y_LN.col(2);
    orig_Y_LN.col(2) = orig_Y_LN.col(1);
    orig_Y_LN.col(1) = temp; 
  }
  

  ///////////////////////////////////////////////////////////////////////////////////////////
  bool bConstVar = (CA->disttype == distribution::normal);
  int  degree    = CA->parms - 2 - (CA->disttype == distribution::normal_ncv );
  double neg_like; 
  // 
  ///////////////////////////////////////////////////////////////////////////////////////////
  normalPOLYNOMIAL_BMD_NC  likelihood_npoly(orig_Y , X, suff_stat, bConstVar, degree);
  normalHILL_BMD_NC  likelihood_nhill(orig_Y , X, suff_stat, bConstVar, 0);
  normalPOWER_BMD_NC likelihood_power(orig_Y , X, suff_stat, bConstVar, 0);
  normalFUNL_BMD_NC  likelihood_funl(orig_Y, X, suff_stat, bConstVar, 0);
  normalEXPONENTIAL_BMD_NC likelihood_nexp5U(orig_Y , X, suff_stat, bConstVar, NORMAL_EXP5_UP);
  normalEXPONENTIAL_BMD_NC likelihood_nexp3U(orig_Y , X, suff_stat, bConstVar, NORMAL_EXP3_UP);
  normalEXPONENTIAL_BMD_NC likelihood_nexp5D(orig_Y , X, suff_stat, bConstVar, NORMAL_EXP5_DOWN);
  normalEXPONENTIAL_BMD_NC likelihood_nexp3D(orig_Y , X, suff_stat, bConstVar, NORMAL_EXP3_DOWN);
  
  lognormalHILL_BMD_NC  likelihood_lnhill(Y_LN, X, suff_stat, 0);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp5U(Y_LN, X, suff_stat, NORMAL_EXP5_UP);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp3U(Y_LN, X, suff_stat, NORMAL_EXP3_UP);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp5D(Y_LN, X, suff_stat, NORMAL_EXP5_DOWN);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp3D(Y_LN, X, suff_stat, NORMAL_EXP3_DOWN);
  
  Eigen::MatrixXd mean; 
  Eigen::MatrixXd var; 
  Eigen::MatrixXd theta(MR->nparms,1); 
  for (int i=0; i < MR->nparms; i++){
    theta(i,0) = MR->parms[i];  
  }
 
  if (CA->disttype == distribution::log_normal){
    switch (CA->model){
      case cont_model::hill:
          mean = likelihood_lnhill.mean(theta,myX); 
          var  = likelihood_lnhill.variance(theta,myX); 
          neg_like = likelihood_lnhill.negLogLikelihood(theta); 
          break; 
      case cont_model::exp_3:
        if (CA->isIncreasing){
          mean = likelihood_lnexp3U.mean(theta,myX); 
          var  = likelihood_lnexp3U.variance(theta,myX);
          neg_like = likelihood_lnexp3U.negLogLikelihood(theta); 
        }else{
          mean = likelihood_lnexp3D.mean(theta,myX); 
          var  = likelihood_lnexp3D.variance(theta,myX);
          neg_like = likelihood_lnexp3D.negLogLikelihood(theta); 
        }
        break; 
      case cont_model::exp_5:
        if (CA->isIncreasing){
          mean = likelihood_lnexp5U.mean(theta,myX); 
          var  = likelihood_lnexp5U.variance(theta,myX);
          neg_like = likelihood_lnexp5U.negLogLikelihood(theta); 
        }else{
          mean = likelihood_lnexp5D.mean(theta,myX); 
          var  = likelihood_lnexp5D.variance(theta,myX);
          neg_like = likelihood_lnexp5D.negLogLikelihood(theta); 
        }
        break; 
      break; 
    }
  }else{

    switch (CA->model){
      case cont_model::funl:
          mean = likelihood_funl.mean(theta,myX); 
          var  = likelihood_funl.variance(theta,myX);
          neg_like = likelihood_funl.negLogLikelihood(theta); 
          break; 
      case cont_model::power:
          mean = likelihood_power.mean(theta,myX); 
          var  = likelihood_power.variance(theta,myX); 
          neg_like = likelihood_power.negLogLikelihood(theta); 
          break; 
      case cont_model::polynomial:
          mean = likelihood_npoly.mean(theta,myX); 
          var  = likelihood_npoly.variance(theta,myX); 
          neg_like = likelihood_npoly.negLogLikelihood(theta); 
          break; 
      case cont_model::hill:
          mean = likelihood_nhill.mean(theta,myX); 
          var  = likelihood_nhill.variance(theta,myX); 
          neg_like = likelihood_nhill.negLogLikelihood(theta); 
          break; 
      case cont_model::exp_3:
        if (CA->isIncreasing){
          mean = likelihood_nexp3U.mean(theta,myX); 
          var  = likelihood_nexp3U.variance(theta,myX);
          neg_like = likelihood_nexp3U.negLogLikelihood(theta); 
        }else{
          mean = likelihood_nexp3D.mean(theta,myX); 
          var  = likelihood_nexp3D.variance(theta,myX);
          neg_like = likelihood_nexp3D.negLogLikelihood(theta); 
        }
        break; 
      case cont_model::exp_5:
        if (CA->isIncreasing){
          mean = likelihood_nexp5U.mean(theta,myX); 
          var  = likelihood_nexp5U.variance(theta,myX);
          neg_like = likelihood_nexp5U.negLogLikelihood(theta); 
        }else{
          mean = likelihood_nexp5D.mean(theta,myX); 
          var  = likelihood_nexp5D.variance(theta,myX);
          neg_like = likelihood_nexp5D.negLogLikelihood(theta); 
        }
        break; 
      break; 
    
    }
    
  }

  
  for (int i = 0; i < expected->n ; i++){
    expected->expected[i] = mean(i,0);
    expected->sd[i] = pow(var(i,0),0.5);
    expected->like = neg_like; 
  }
  
}
