#include <analysis_of_deviance.h>
#include <mcmc_analysis.h>
#include <bmd_calculate.h>
#include <IDPriorMCMC.h>
 
#include <chrono>
#include <algorithm>
#include <vector>
#include <limits>
#include <set>
#include <numeric>

void inverse_transform_dose(continuous_model_result *model){
  if (model){
    model->bmd = sinh(model->bmd); 
    for (int i = 0; i< model->dist_numE; i ++){
      double temp = double(i)/double(model->dist_numE); 
      model->bmd_dist[i] = sinh(model->bmd_dist[i]);     // BMD @ probability
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

// Eigen::MatrixXd init_funl_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
//   std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
//   std::sort(vec.begin(), vec.end());
//   vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
//   //udoses = vec; // this should be the unique dose group
  
//   Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
//   prior(0,1) = betas(0,0); 
//   double max_d = vec[vec.size()-1]; 
//   double max_r = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d);
//   prior(1,1)   = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d - prior(0,1))/max_d; 
//   prior(2,1)   = max_r; 
//   prior(3,1)   = 0.5;
//   prior(4,1)   = 1;
//   prior(5,1)   = 0.75;
//   prior(6,1)   = 1;
  
//   for (int i = 0; i < 7; i++){
//     if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3); 
//     if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
//   }
  
  
//   return prior; 
// }

// Eigen::MatrixXd init_test4_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){

//   double minDose = X.minCoeff();
//   double maxDose = X.maxCoeff();
//   double init = 0;
//   int nmin = 0, nmax = 0;

//   for (int i = 0; i < X.rows(); i++){
//     if (X(i,0)==minDose){
//       nmin++;
//       init += Y_N(i,0);
//     }
//   }
//   init *= init/double(nmin);
//   prior(0,1) = init;

//   init = 0;
//   for (int i = 0; i < X.rows(); i++){
//     if (X(i,0)==maxDose){
//       nmax++;
//       init += Y_N(i,0);
//     }
//   }
//   init *= init/double(nmax);

//   prior(2,1)   =  init / prior(0,1);
//   prior(1,1)   = 0.0001*maxDose;
//   prior(3,1)   = 5;

//   //make sure the starting point is within bounds; if not, put on boundary
//   for(int i = 0; i < 4; i++){
// 	  if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3);
// 	  if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
//   }
//   //cerr << prior << endl;
//   return prior;
// }

// Eigen::MatrixXd init_test4_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
//   Y_LN.col(0) = exp(Y_LN.col(0).array());
//   if (Y_LN.cols() ==3 ){
//     Y_LN.col(1) = exp(Y_LN.col(1).array());
//   }
//   return init_test4_nor(Y_LN,  X, prior);

// }

// Eigen::MatrixXd init_test5_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){

//   double minDose = X.minCoeff();
//   double maxDose = X.maxCoeff();
//   double init = 0;
//   int nmin = 0, nmax = 0;

//   for (int i = 0; i < X.rows(); i++){
//     if (X(i,0)==minDose){
//       nmin++;
//       init += Y_N(i,0);
//     }
//   }
//   init *= init/double(nmin);
//   prior(0,1) = init;

//   init = 0;
//   for (int i = 0; i < X.rows(); i++){
//     if (X(i,0)==maxDose){
//       nmax++;
//       init += Y_N(i,0);
//     }
//   }
//   init *= init/double(nmax);

//   prior(2,1)   =  init / prior(0,1);
//   prior(1,1)   = 0.0001*maxDose;
//   prior(3,1)   = 5;
//   prior(4,1)   = 2;

//   //make sure the starting point is within bounds; if not, put on boundary
//   for(int i = 0; i < 5; i++){
// 	  if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3);
// 	  if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
//   }
//   //cerr << prior << endl;
//   return prior;
// }

// Eigen::MatrixXd init_test5_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
//   Y_LN.col(0) = exp(Y_LN.col(0).array());
//   if (Y_LN.cols() ==3 ){
//     Y_LN.col(1) = exp(Y_LN.col(1).array());
//   }
//   return init_test5_nor(Y_LN,  X, prior);

// }


// Eigen::MatrixXd init_hill_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
//   std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
//   std::sort(vec.begin(), vec.end());
//   vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
//   //udoses = vec; // this should be the unique dose group
//   double minDose = X.minCoeff(); 
//   double maxDose = X.maxCoeff(); 
//   double init = 0; 
//   int nmin = 0, nmax = 0;  
  
//   for (int i = 0; i < X.rows(); i++){
//     if (X(i,0)==minDose){
//       nmin++; 
//       init += Y_N(i,0); 
//     }
//   }
//   init *= init/double(nmin); 
  
//   Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
//   prior(0,1) = init; 
//   init = 0;
//   for (int i = 0; i < X.rows(); i++){
//     if (X(i,0)==maxDose){
//       nmax++; 
//       init += Y_N(i,0); 
//     }
//   }
//   init *= init/double(nmin); 
  
//   prior(1,1)   =  (init - prior(0,1))/(maxDose-minDose); 
//   prior(2,1)   = 0;//0.0001*maxDose; 
//   prior(3,1)   = 10;
  
//   if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
//   if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
//   if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
//   if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
//   //cerr << prior << endl; 
//   return prior; 
// }


// Eigen::MatrixXd init_pow_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
//   std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
//   std::sort(vec.begin(), vec.end());
//   vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
//   //udoses = vec; // this should be the unique dose group
  
//   Eigen::MatrixXd betas = powerSearchRegression(Y_N,  X);
//   prior(0,1)   = betas(0,0); 
//   prior(1,1)   = betas(1,0);  
//   prior(2,1)   = betas(2,0);
  
//   if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
//   if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
//   if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
//   if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  
//   if (prior(2,1) < prior(1,3)) prior(2,1) = prior(1,3); 
//   if (prior(2,1) > prior(1,4)) prior(2,1) = prior(1,4);
  
//   return prior; 
// }

// Eigen::MatrixXd init_hill_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
//   Y_LN.col(0) = exp(Y_LN.col(0).array());
//   if (Y_LN.cols() ==3 ){
//     Y_LN.col(1) = exp(Y_LN.col(1).array());
//   }
//   return init_hill_nor(Y_LN,  X, prior); 
  
// }


// Eigen::MatrixXd init_exp_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
//   std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
//   std::sort(vec.begin(), vec.end());
//   vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
//   //udoses = vec; // this should be the unique dose group
  
//   Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
//   prior(0,1) = betas(0,0); 
//   double max_d = vec[vec.size()-1]; 
//   double max_r = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d);
//   prior(2,1)   = log(0.001);  
//   double temp = max_r/prior(0,1);
  
//   temp =  -(temp-exp(prior(2,1)))/(exp(prior(2,1))-1.0);
  
//   prior(1,1) = 0.05; 
//   prior(3,1)   = 2.5; 
  
//   if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
//   if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
//   if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
//   if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  
  
//   return prior; 
// }

// Eigen::MatrixXd init_exp_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
//  //right here
//   Y_LN.col(0) = exp(Y_LN.col(0).array());
//   if (Y_LN.cols() ==3 ){
//       Y_LN.col(1) = exp(Y_LN.col(1).array());
//   }
//   return init_exp_nor(Y_LN, X, prior); 
// }

// Eigen::MatrixXd init_poly(Eigen::MatrixXd Y, Eigen::MatrixXd tX, 
//                           Eigen::MatrixXd prior, int deg = 2){

//   Eigen::MatrixXd X = Eigen::MatrixXd::Ones(tX.rows(),deg+1);
//   Eigen::MatrixXd W = Eigen::MatrixXd::Identity(tX.rows(),tX.rows());
 
//   for (int i = 0; i < X.rows(); i++){
//     if (Y.cols()>1){
//       W(i,i) = Y(i,2)/Y(i,1)*Y(i,1); // Weights: \sigma^2/N
//     }
//     for (int j = 1; j < X.cols(); j++){
//       X(i,j) = pow(tX(i,0),j);  
//     }
//   }
//   Eigen::MatrixXd B = Eigen::MatrixXd::Ones(deg+1,1);
//   B = X.transpose()*W*X;
//   B = B.inverse()*X.transpose()*W*Y.col(0); 
//   for(int i = 0; i < B.rows(); i++){
//     if ( B(i,0) < prior(i,3) ){
//       prior(i,1) = prior(i,3); 
//     } else if (B(i,0) > prior(i,4)){
//       prior(i,1) = prior(i,4); 
//     }else{
//       prior(i,1) = B(i,0);
//     }
//   } 
  
//   return prior; 
// }

// /*initialize_mle
//  * This function is for MLE optimization it takes the data/model type and then tries 
//  * to start the initializer on reasonable initial values. These values will then be fed
//  * to the optimizer.
//  * OUTPUT: new prior vector with the new initial values put in the mean column. 
//  */
// Eigen::MatrixXd initialize_model(Eigen::MatrixXd Y_N, Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, 
//                                  Eigen::MatrixXd prior, distribution data_dist, cont_model model){
  
//   Eigen::MatrixXd retVal = prior; 
//   int deg = 0; 
//   switch (model){
//   case cont_model::funl:
//     retVal = init_funl_nor(Y_N, X, prior);
//     break; 
//   case cont_model::hill:
//     retVal = distribution::log_normal==data_dist ? init_hill_lognor(Y_LN, X, prior):
//     init_hill_nor(Y_N, X, prior); 
//     break; 
//   case cont_model::exp_aerts: case cont_model::invexp_aerts: case cont_model::hill_aerts: case cont_model::lognormal_aerts:
//   case cont_model::logistic_aerts: case cont_model::probit_aerts: case cont_model::LMS: case cont_model::gamma_efsa:
//     retVal = distribution::log_normal==data_dist ? init_test4_lognor(Y_LN, X, prior):
//     init_test4_nor(Y_N, X, prior);
//     break;
//   case cont_model::gamma_aerts: case cont_model::invgamma_aerts: case cont_model::lomax_aerts:
//   case cont_model::invlomax_aerts: case cont_model::logskew_aerts: case cont_model::invlogskew_aerts:
// 	retVal = distribution::log_normal==data_dist ? init_test5_lognor(Y_LN, X, prior):
// 	init_test5_nor(Y_N, X, prior);
// 	break;
//   case cont_model::exp_3:
//   case cont_model::exp_5:
//     retVal = distribution::log_normal==data_dist ? init_exp_lognor(Y_LN, X, prior):
//                                                    init_exp_nor(Y_N, X, prior); 
//     break; 
//   case cont_model::power: 

//     retVal = init_pow_nor( Y_N,  X,  prior);
//     break; 
//   case cont_model::polynomial:
//     /*initialize at Least Squares inv(X.t()*X)*X.t()*Y)
//      * 
//      */

//     deg = distribution::normal_ncv == data_dist? prior.rows() - 3: prior.rows() - 2;   
//     retVal = distribution::log_normal==data_dist ? init_poly(Y_LN, X, prior,deg):
//                                                    init_poly(Y_N, X, prior,deg); 
//     break; 
//   default: 
//     // this is a different model that shouldn't have MLE fits
//     // so we don't do anything
//     break; 
//   }

//   return retVal.col(1);  
// }


// double compute_lognormal_dof(Eigen::MatrixXd Y,Eigen::MatrixXd X, Eigen::MatrixXd estimate, 
//                              bool is_increasing, bool suff_stat, Eigen::MatrixXd prior, 
//                              cont_model CM){
//   double DOF = 0; 
//   Eigen::MatrixXd Xd; 
//   Eigen::MatrixXd cv_t; 
//   Eigen::MatrixXd pr; 
//   Eigen::MatrixXd temp(X.rows(),3);
//   Eigen::MatrixXd subBlock(3,3); 
//   Eigen::MatrixXd temp_estimate(estimate.rows() + 1,1); 
  
//   switch(CM){
//   case cont_model::hill:
//     Xd = X_gradient_cont<lognormalHILL_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4); 
//     cv_t = X_cov_cont<lognormalHILL_BMD_NC>(estimate,Y,X,suff_stat); 
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     pr = pr.block(0,0,4,4); 
    
//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0; 
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr; 
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
    
//     break; 
//   case cont_model::exp_aerts:
//     Xd = X_gradient_cont<lognormalEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::LMS:
//     Xd = X_gradient_cont<lognormalLMS_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalLMS_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::gamma_efsa:
//     Xd = X_gradient_cont<lognormalGAMMA_efsa_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalGAMMA_efsa_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invexp_aerts:
//     Xd = X_gradient_cont<lognormalIEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalIEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::gamma_aerts:
//     Xd = X_gradient_cont<lognormalGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),5);
//     cv_t = X_cov_cont<lognormalGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,5,5);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 5.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invgamma_aerts:
//     Xd = X_gradient_cont<lognormalIGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),5);
//     cv_t = X_cov_cont<lognormalIGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,5,5);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 5.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::hill_aerts:
//     Xd = X_gradient_cont<lognormalHILL_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalHILL_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::lomax_aerts:
//     Xd = X_gradient_cont<lognormalLOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),5);
//     cv_t = X_cov_cont<lognormalLOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,5,5);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 5.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invlomax_aerts:
//     Xd = X_gradient_cont<lognormalILOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),5);
//     cv_t = X_cov_cont<lognormalILOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,5,5);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 5.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::lognormal_aerts:
//     Xd = X_gradient_cont<lognormalLOGNORMAL_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalLOGNORMAL_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::logskew_aerts:
//     Xd = X_gradient_cont<lognormalLOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),5);
//     cv_t = X_cov_cont<lognormalLOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,5,5);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 5.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invlogskew_aerts:
//     Xd = X_gradient_cont<lognormalILOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),5);
//     cv_t = X_cov_cont<lognormalILOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,5,5);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 5.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::logistic_aerts:
//     Xd = X_gradient_cont<lognormalLOGISTIC_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalLOGISTIC_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::probit_aerts:
//     Xd = X_gradient_cont<lognormalPROBIT_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     Xd = Xd.block(0,0,Xd.rows(),4);
//     cv_t = X_cov_cont<lognormalPROBIT_aerts_BMD_NC>(estimate,Y,X,suff_stat);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     pr = pr.block(0,0,4,4);

//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0;
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::exp_3:
    
//     temp_estimate << estimate(0,0) , estimate(1,0) , 1.0 , estimate.block(2,0,estimate.rows()-2,1); 
//     if (is_increasing){
//       Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_UP);
//       cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_UP);
//     }else{
//       Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_DOWN);
//       cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_DOWN);
//     }
    
//     temp << Xd.col(0) , Xd.col(1), Xd.col(3);  
//     Xd = temp; 
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     subBlock << pr(0,0), pr(0,1), pr(0,3),
//                 pr(1,0), pr(1,1), pr(1,3),
//                 pr(3,0), pr(3,1), pr(3,3);
    
//     if( fabs(subBlock.diagonal().array().sum()) ==0){
//       DOF = 3; 
//     } else{
//       pr   = Xd.transpose()*cv_t*Xd + subBlock; 
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
//     break;
//   case cont_model::exp_5: 
//   default:
//     if (is_increasing){
//       Xd = X_gradient_cont<lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,NORMAL_EXP5_UP);
//       cv_t = X_cov_cont< lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat, NORMAL_EXP5_UP);
//     }else{
//       Xd = X_gradient_cont<lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,NORMAL_EXP5_DOWN);
//       cv_t = X_cov_cont< lognormalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat, NORMAL_EXP5_DOWN);
//     }
  
//     Eigen::MatrixXd temp_Xd = Xd.block(0,0,Xd.rows(),4); 
//     Xd = temp_Xd; 
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     temp_Xd  = pr.block(0,0,4,4); 
//     pr = temp_Xd; 
//     if( fabs(pr.diagonal().array().sum()) == 0){
//       DOF = 4.0; 
//     }else{
//       pr   = Xd.transpose()*cv_t*Xd + pr; 
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
//     break; 
//   }

//   return DOF; 
// }

// double compute_normal_dof(Eigen::MatrixXd Y,Eigen::MatrixXd X, Eigen::MatrixXd estimate, 
//                           bool is_increasing, bool suff_stat,  bool CV, Eigen::MatrixXd prior,
//                           cont_model CM,int degree){
//   double DOF = 0; 
//   Eigen::MatrixXd Xd; 
//   Eigen::MatrixXd cv_t; 
//   Eigen::MatrixXd pr; 
//   Eigen::MatrixXd temp(X.rows(),3);
//   Eigen::MatrixXd subBlock(3,3); 
 
//   int offset = CV? 1:2; 
//   Eigen::MatrixXd temp_estimate(estimate.rows() + 1,1); 
//   Eigen::MatrixXd temp_block(1,1); 

//   switch(CM){
//   case cont_model::polynomial:
    
//     Xd = X_gradient_cont_norm<normalPOLYNOMIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,degree);
//     temp_block = Xd.block(0,0,Xd.rows(),estimate.rows() - offset); 
//     Xd = temp_block; 
//     cv_t = X_cov_cont_norm<normalPOLYNOMIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,degree); 
    
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
    
//     temp_block = pr.block(0,0,estimate.rows() - offset,estimate.rows() - offset); 
//     pr = temp_block; 
    
//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = pr.diagonal().size(); 
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr; 
//       pr = temp_block; 
//       temp_block = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       Xd = temp_block; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
//     break; 
//   case cont_model::hill:
//     Xd = X_gradient_cont_norm<normalHILL_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4); 
//     Xd = temp_block; 
//     cv_t = X_cov_cont_norm<normalHILL_BMD_NC>(estimate,Y,X,suff_stat,CV); 
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     temp_block  =    pr.block(0,0,4,4); 
//     pr = temp_block; 
    
//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0; 
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr; 
//       pr = temp_block; 
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
    
//     break; 
//   case cont_model::exp_aerts:
//     Xd = X_gradient_cont_norm<normalEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::LMS:
//     Xd = X_gradient_cont_norm<normalLMS_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalLMS_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::gamma_efsa:
//     Xd = X_gradient_cont_norm<normalGAMMA_efsa_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalGAMMA_efsa_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invexp_aerts:
//     Xd = X_gradient_cont_norm<normalIEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalIEXP_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::hill_aerts:
//     Xd = X_gradient_cont_norm<normalHILL_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalHILL_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::lognormal_aerts:
//     Xd = X_gradient_cont_norm<normalLOGNORMAL_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalLOGNORMAL_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::logistic_aerts:
//     Xd = X_gradient_cont_norm<normalLOGISTIC_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalLOGISTIC_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::probit_aerts:
//     Xd = X_gradient_cont_norm<normalPROBIT_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),4);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalPROBIT_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,4,4);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 4.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::gamma_aerts:
//     Xd = X_gradient_cont_norm<normalGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),5);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,5,5);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 5.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invgamma_aerts:
//     Xd = X_gradient_cont_norm<normalIGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),5);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalIGAMMA_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,5,5);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 5.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::lomax_aerts:
//     Xd = X_gradient_cont_norm<normalLOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),5);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalLOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,5,5);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 5.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invlomax_aerts:
//     Xd = X_gradient_cont_norm<normalILOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),5);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalILOMAX_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,5,5);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 5.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::logskew_aerts:
//     Xd = X_gradient_cont_norm<normalLOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),5);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalLOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,5,5);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 5.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::invlogskew_aerts:
//     Xd = X_gradient_cont_norm<normalILOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     temp_block  = Xd.block(0,0,Xd.rows(),5);
//     Xd = temp_block;
//     cv_t = X_cov_cont_norm<normalILOGSKEW_aerts_BMD_NC>(estimate,Y,X,suff_stat,CV);
//     pr   =  X_logPrior<IDPrior>(estimate,prior);
//     temp_block  =    pr.block(0,0,5,5);
//     pr = temp_block;

//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF = 5.0;
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr;
//       pr = temp_block;
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t;
//       DOF =  Xd.diagonal().array().sum();
//     }

//     break;
//   case cont_model::exp_3:
    
//     temp_estimate << estimate(0,0) , estimate(1,0) , 1.0 , estimate.block(2,0,estimate.rows()-2,1); 
//     if (is_increasing){
//       Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_UP);
//       cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_UP);
//     }else{
//       Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat,NORMAL_EXP3_DOWN);
//       cv_t = X_cov_cont_norm<normalEXPONENTIAL_BMD_NC>(temp_estimate,Y,X,suff_stat, NORMAL_EXP3_DOWN);
//     }

//     temp << Xd.col(0) , Xd.col(1), Xd.col(3);  
//     Xd = temp; 
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     subBlock << pr(0,0), pr(0,1), pr(0,3),
//                 pr(1,0), pr(1,1), pr(1,3),
//                 pr(3,0), pr(3,1), pr(3,3);
    
//     if( fabs(subBlock.diagonal().array().sum()) ==0){
//       DOF = 3; 
//     } else{
//       pr   = Xd.transpose()*cv_t*Xd + subBlock; 
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
//     break; 
//   case cont_model::exp_5: 
//     if (is_increasing){
//       Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_UP);
//       cv_t = X_cov_cont_norm< normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_UP);
//     }else{
//       Xd = X_gradient_cont_norm<normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_DOWN);
//       cv_t = X_cov_cont_norm< normalEXPONENTIAL_BMD_NC>(estimate,Y,X,suff_stat,CV,NORMAL_EXP5_DOWN);
//     }
    
//     temp_block = Xd.block(0,0,Xd.rows(),3); 
//     Xd = temp_block; 
    
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     temp_block = pr.block(0,0,3,3); 
//     pr = temp_block; 
//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF =4.0; 
//     } else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr; 
//       pr = temp_block; 
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
    
//     break; 
//   case cont_model::power: 
//   default:
//     Xd = X_gradient_cont_norm<normalPOWER_BMD_NC>(estimate,Y,X,CV,suff_stat);
//     cv_t = X_cov_cont_norm<normalPOWER_BMD_NC>(estimate,Y,X,CV,suff_stat);
    
//     temp_block = Xd.block(0,0,Xd.rows(),3); 
//     Xd = temp_block; 
//     pr   =  X_logPrior<IDPrior>(estimate,prior); 
//     temp_block = pr.block(0,0,3,3); 
//     pr = temp_block; 
    
//     if( fabs(pr.diagonal().array().sum()) ==0){
//       DOF =3.0; 
//     }else{
//       temp_block   = Xd.transpose()*cv_t*Xd + pr; 
//       pr = temp_block;  
//       Xd = Xd*pr.inverse()*Xd.transpose()*cv_t; 
//       DOF =  Xd.diagonal().array().sum(); 
//     }
    
//     break;
//   }   
  
//   return DOF + offset; 
  
// }

// mcmcSamples mcmc_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
//                             Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
//                             bool is_increasing, 
//                             double bmrf,   double bk_prob, 
//                             double alpha,  int samples, int burnin,  
//                             Eigen::MatrixXd initV) {
   
//   bool suff_stat = Y.cols() == 1? false:true; 
//   std::vector<bool> fixedB(prior.rows());
//   std::vector<double> fixedV(prior.rows());
//   for (int i = 0; i < prior.rows(); i++) {
//     fixedB[i] = false;
//     fixedV[i] = 0.0;
//   }

//   mcmcSamples a;
//   int adverseR; 
//   switch (CM)
//   {
//   case cont_model::hill:
// #ifdef R_COMPILATION 
//     //cout << "Running Hill Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalHILL_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV); 
//     break; 
//   case cont_model::exp_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Exponential (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXP_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::LMS:
// #ifdef R_COMPILATION
//     //cout << "Running LMS Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalLMS_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::gamma_efsa:
// #ifdef R_COMPILATION
//     //cout << "Running Gamma-EFSA Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalGAMMA_efsa_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::invexp_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Inverse Exponential (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalIEXP_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::gamma_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Gamma (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalGAMMA_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::invgamma_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Inverse Gamma (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalIGAMMA_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::hill_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Hill (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalHILL_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::lomax_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Lomax (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalLOMAX_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::invlomax_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Inverse Lomax (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalILOMAX_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::lognormal_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Lognormal (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalLOGNORMAL_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::logskew_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Logskew (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalLOGSKEW_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::invlogskew_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Inverse Logskew (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalILOGSKEW_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::logistic_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Logistic (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalLOGISTIC_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::probit_aerts:
// #ifdef R_COMPILATION
//     //cout << "Running Probit (Aerts)  Model Log-Normality Assumption using MCMC." << endl;
// #endif
//      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalPROBIT_aerts_BMD_NC, IDPriorMCMC>
//                                     (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                      bk_prob,suff_stat,bmrf, riskType, alpha,
//                                      samples,0,initV);
//     break;
//   case cont_model::exp_3:
//       adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
// #ifdef R_COMPILATION 
//       //cout << "Running Exponential 3 Model Log-Normality Assumption using MCMC." << endl;
// #endif
//       a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXPONENTIAL_BMD_NC, IDPriorMCMC>
//                                               (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                               bk_prob,suff_stat,bmrf, riskType, alpha,
//                                               samples,adverseR,initV);
//       // remove the third entry
//       removeRow(a.map_cov, 2);
//       removeCol(a.map_cov, 2);
//       removeRow(a.map_estimate, 2);
//     break; 
//   case cont_model::exp_5:
//   default: 
//       adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
// #ifdef R_COMPILATION 
//       //cout << "Running Exponential 5 Model Log-Normality Assumption using MCMC." << endl;
// #endif
//       a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXPONENTIAL_BMD_NC, IDPriorMCMC>
//                                               (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                bk_prob,suff_stat,bmrf, riskType, alpha,
//                                                samples,adverseR,initV);
//     break; 
    
//   }
  
//   return a; 
// }

// void inverse_transform_dose(bmd_analysis_MCMC *b){
  
//   if (b){
//     for(unsigned int i= 0; i < b->samples;  i++){
//       b->BMDS[i] = sinh(b->BMDS[i]); 
//     }
//   }
  
// }

// void estimate_normal_variance(continuous_analysis *CA,
//                               double *v_c, double *v_nc, double *v_pow){
  
//   // standardize the data
//   int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
//   Eigen::MatrixXd Y(n_rows,n_cols); 
//   Eigen::MatrixXd X(n_rows,1); 
//   // copy the origional data
  
//   for (int i = 0; i < n_rows; i++){
//     Y(i,0) = CA->Y[i]; 
//     X(i,0) = CA->doses[i]; 
//     if(CA->suff_stat){
//       Y(i,2) = CA->sd[i]; 
//       Y(i,1) = CA->n_group[i]; 
//     }
//   }
  
//   double divisor = get_divisor( Y,  X); 
//   double  max_dose = X.maxCoeff(); 
  
//   Eigen::MatrixXd orig_Y = Y, orig_Y_LN = Y; 
//   Eigen::MatrixXd orig_X = X; 
  
//   Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
//   Eigen::MatrixXd Y_LN, Y_N;
//   bool can_be_suff = true; 
//   if (Y.cols() == 1){ 
//     //individual data
//     can_be_suff = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
//   }else{
    
//     SSTAT     = cleanSuffStat(Y,X,false,false); 
//     SSTAT_LN  = cleanSuffStat(Y,X,true,false);
//     UX = X; 
//   }

  
//   if(!can_be_suff){
//     *v_c   =  std::numeric_limits<double>::infinity();
//     *v_nc  =  std::numeric_limits<double>::infinity();
//     *v_pow =  std::numeric_limits<double>::infinity();
//     return;   
//   }else{
//      Y_N = SSTAT; 
//      Eigen::MatrixXd temp = Y_N.col(2);
//      Y_N.col(2) = Y_N.col(1);
//      Y_N.col(1) = temp; 
//      variance_fits(Y_N, UX, 
//                     can_be_suff, 
//                     v_c, v_nc, v_pow);
    
//     return; 
//   }
// }





// mcmcSamples mcmc_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
//                         Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
//                         bool is_increasing, bool bConstVar,
//                         double bmrf,   double bk_prob, 
//                         double alpha, int samples,
//                         int burnin, Eigen::MatrixXd initV,
//                         int degree = 2); 
// mcmcSamples mcmc_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
//                          Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
//                          bool is_increasing, bool bConstVar,
//                          double bmrf,   double bk_prob, 
//                          double alpha, int samples,
//                          int burnin, Eigen::MatrixXd initV,
//                          int degree) {
  
//   bool suff_stat = Y.cols() == 1? false:true; 
//   std::vector<bool> fixedB(prior.rows());
//   std::vector<double> fixedV(prior.rows());
  
//   for (int i = 0; i < prior.rows(); i++) {
//     fixedB[i] = false;
//     fixedV[i] = 0.0;
//   }
  
//   mcmcSamples a;
//   int adverseR; 
//   switch (CM)
//   {
//   case cont_model::funl:
// #ifdef R_COMPILATION 
//     if (bConstVar){
//       //cout << "Running FUNL Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running FUNL Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif  
//     a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalFUNL_BMD_NC, IDPriorMCMC>
//                                             (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                              bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                              samples,0,initV); 
//     break;   
    
//   case cont_model::hill:
// #ifdef R_COMPILATION 
//     if (bConstVar){
//       //cout << "Running Hill Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Hill Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif
    
//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalHILL_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV); 
//     break; 
//   case cont_model::exp_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Exponential (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Exponential (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXP_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::LMS:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running LMS Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running LMS Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalLMS_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::gamma_efsa:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Gamma-EFSA Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Gamma-EFSA Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalGAMMA_efsa_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::invexp_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Inverse Exponential (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Inverse Exponential (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalIEXP_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::gamma_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Gamma (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Gamma (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalGAMMA_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::invgamma_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Inverse Gamma (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Inverse Gamma (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalIGAMMA_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::hill_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Hill (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Hill (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalHILL_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::lomax_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Lomax (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Lomax (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalLOMAX_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::invlomax_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Inverse Lomax (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Inverse Lomax (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalILOMAX_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::lognormal_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Log-normal (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Log-normal (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalLOGNORMAL_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::logskew_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Log-skew-normal (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Log-skew-normal (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalLOGSKEW_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::invlogskew_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Inverse Log-skew-normal (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Inverse Log-skew-normal (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalILOGSKEW_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::logistic_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Logistic (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Logistic (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalLOGISTIC_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::probit_aerts:
// #ifdef R_COMPILATION
//     if (bConstVar){
//       //cout << "Running Probit (Aerts)  Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Probit (Aerts)  Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif

//       a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPROBIT_aerts_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break;
//   case cont_model::exp_3:
//     adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
// #ifdef R_COMPILATION 
//     if (bConstVar){
//       //cout << "Running Exponential 3 Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Exponential 3 Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif
//     a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,adverseR,initV);
//     // remove the third entry
//     removeRow(a.map_cov, 2);
//     removeCol(a.map_cov, 2);
//     removeRow(a.map_estimate, 2);
    
//     break; 
//   case cont_model::exp_5:
//     adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
// #ifdef R_COMPILATION 
//     if (bConstVar){
//       //cout << "Running Exponential 5 Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Exponential 5 Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif
//     a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDPriorMCMC>
//                                                 (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                                  bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                                  samples,0,initV);
//     break; 
//   case cont_model::power:
 
// #ifdef R_COMPILATION 
//     if (bConstVar){
//       //cout << "Running Power Model Normality Assumption using MCMC." << endl;
//     }else{
//       //cout << "Running Power Model Normality-NCV Assumption using MCMC." << endl;
//     }
// #endif
//     a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPOWER_BMD_NC, IDPriorMCMC>
//                                               (Y,  X, prior, fixedB, fixedV, is_increasing,
//                                               bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//                                               samples,adverseR,initV);
    
  
//   break; 
//   case cont_model::polynomial:
// #ifdef R_COMPILATION 
//   if (bConstVar){
//     //cout << "Running Polynomial Model Normality Assumption using MCMC." << endl;
//   }else{
//     //cout << "Running Polynomial Model Normality-NCV Assumption using MCMC." << endl;
//   }
// #endif

//     a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPOLYNOMIAL_BMD_NC, IDPriorMCMC>
//           (Y,  X, prior, fixedB, fixedV, is_increasing,
//            bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
//            samples,degree,initV);
    
//   default: 
//     break; 
//   }
//   //convert a stuff
//   //
//   return a; 
// }
