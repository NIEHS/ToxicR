/*
 
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
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]

#ifdef R_COMPILATION
//necessary things to run in R    
  #include <RcppEigen.h>
  #include <RcppGSL.h>
using namespace Rcpp;
#else 
  #include <Eigen/Dense>

#endif

#include <gsl/gsl_randist.h>


#include "bmdStruct.h"
#include "continuous_clean_aux.h"
#include <iostream>
#include <numeric> 


using namespace std;


// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]


/*data clean procedures*/ 
 
 std::vector<double> unique_list(Eigen::MatrixXd X){
   
   std::vector<double> uniqueX; 
   
   for (int i = 0; i < X.rows(); i++){
     bool isUnique = true; 
     for (int j = 0; j < uniqueX.size(); j++){
       if (X(i,0)==uniqueX[j]){
         isUnique = false; 
         break; 
       }
     }
     if (isUnique){
       uniqueX.push_back(X(i,0)); 
     }
   }
   
   return uniqueX; 
 }
 
 
Eigen::MatrixXd cleanSuffStat(Eigen::MatrixXd Y, Eigen::MatrixXd X, bool is_logNormal, bool use_divisor){
  double minDose = X.minCoeff(); 
  double divisor = 0; 
  int nmin = 0; 
  
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++; 
      divisor += Y(i,0); 
    }
  }
  
  if (use_divisor){
     divisor = get_divisor(Y,X);  //average background dose
  }else{
     divisor = 1.0; 
  } 
  
  
  
  Y.col(0).array() = Y.col(0).array()/divisor; //divide mean
  Y.col(2).array() = Y.col(2).array()/divisor; //divide sd; 
  if (is_logNormal){
    
    Eigen::MatrixXd t1 = sqrt(log(1.0+pow(Y.col(2).array(),2)/Y.col(0).array())); 
    Eigen::MatrixXd t2 = log(Y.col(0).array())-pow(t1.array(),2)*0.5; 
    Y.col(0) = t2.array(); 
    Y.col(2) = t1.array(); 
  }
  return Y; 
  
}

double get_divisor(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  // find the average of the lowest dose (e.g. background)
  double minDose = X.minCoeff(); 
  double divisor = 0; 
  int nmin = 0; 
  
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0) == minDose){
      nmin++; 
      divisor += Y(i,0); 
    }
  }
  divisor = divisor/double(nmin); 
  
  return  fabs(divisor) < 1 ? 1: fabs(divisor); // return the absolute value of the divisor so we don't 
                        // flip the sign of the dose-response curve. 
}

Eigen::MatrixXd createSuffStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                               bool is_logNormal){
  // get unique element
  std::vector<double> uniqueX = unique_list(X); 
  
  // find the average of the lowest dose (e.g. background)
  double minDose = X.minCoeff(); 
  double divisor = 0; 
  int nmin = 0; 
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++; 
      divisor += Y(i,0); 
    }
  }
  divisor =  get_divisor(X,Y) ;
  
  // build the sufficient statistics
  Eigen::MatrixXd SSTAT(uniqueX.size(),3); 
  for (int i = 0; i < uniqueX.size(); i++){
    std::vector<double> uVals; 
    for (int j = 0; j < Y.rows(); j++){
      if (X(j,0)==uniqueX[i]){
        if (is_logNormal){
          uVals.push_back(log(Y(j,0))); //geometric mean
        }else{
          uVals.push_back(Y(j,0));  // aritmetic mean
        }
      }
    }
    SSTAT(i,0) = accumulate( uVals.begin(), uVals.end(), 0.0) / uVals.size();
    SSTAT(i,1) = uVals.size();
    double sqSum = std::inner_product(uVals.begin(), uVals.end(), uVals.begin(), 0.0);
    SSTAT(i,2) =  sqSum / uVals.size() - SSTAT(i,0)*SSTAT(i,0);
    SSTAT(i,2) =  SSTAT(i,2) * (double(uVals.size())/(double(uVals.size())-1.0)); //n-1 instead of n
    SSTAT(i,2) =  pow(SSTAT(i,2),0.5); 
  }
  
  return SSTAT; 
}
 
// FIXME: CHECK IF WE ARE RESCALING THE VARIANCE PARAMETERS CORRECTLY
Eigen::MatrixXd rescale_parms(Eigen::MatrixXd parms, cont_model model, double max_dose, double bkground,
                              bool is_logNormal, int degree)
  {
 
    switch(model){
      case cont_model::hill:
        parms(0,0) *= bkground; parms(1,0) *= bkground; parms(2,0)*=max_dose; 
        if (!is_logNormal){
          if (parms.rows()==5){
            parms(4,0) += 2*log(bkground); 
          }else{
            parms(5,0) += 2*log(bkground); 
          }
        }
        break; 
      case cont_model::exp_3:
     
        parms(0,0) *= bkground; parms(1,0) *= 1/max_dose; 
        if (!is_logNormal){
          if (parms.rows()== 4){
            parms(3,0) += 2*log(bkground); 
          }else{
            parms(4,0) += 2*log(bkground); 
          }
        }
       
        break; 
      case cont_model::exp_5:
        
        parms(0,0) *= bkground; parms(1,0) *= 1/max_dose; 
        if (!is_logNormal){
            if (parms.rows()==5){
              parms(4,0) += 2*log(bkground); 
            }else{
              parms(5,0) += 2*log(bkground); 
            }
        }
        break; 
        
      case cont_model::power: 
        parms(0,0) *= bkground; parms(1,0) *= bkground; 
        parms(1,0) *= pow(1/max_dose,parms(2,0)); 
        if (!is_logNormal){
          if (parms.rows()==4){
            parms(3,0) += 2*log(bkground); 
          }else{
            parms(4,0) += 2*log(bkground); 
          }
        }
        break; 
      case cont_model::polynomial:
      
      for (int i = 1; i <= degree; i++){
        parms(i,0) *= pow(1/max_dose,i); 
      }
      if (!is_logNormal){
        parms(parms.rows()-1,0) += 2*log(bkground); 
      }
     
      break; 
    case cont_model::funl:
      parms(0,0) *= bkground; parms(1,0) *= bkground; 
      parms(2,0) *= max_dose; 
      parms(3,0) *= max_dose; 
      parms(4,0) *= max_dose; 
      parms(5,0) += 2*log(abs(1/max_dose));   
      if (!is_logNormal){
        if (parms.rows()==7){
          parms(6,0) += 2*log(bkground); 
        }else{
          parms(7,0) += 2*log(bkground); 
        }
      }
      break; 
      default:
        // don't know what to do so we return
        break; 
    }
    
    return parms; 
    
  }
 
 
 Eigen::MatrixXd rescale_cov_matrix(Eigen::MatrixXd COV, 
									Eigen::MatrixXd parms, cont_model model,
									double max_dose, double bkground,
									bool is_logNormal,int degree)
  {
   
    Eigen::MatrixXd scaleMatrix = Eigen::MatrixXd::Identity(COV.rows(), COV.cols());
    switch(model){
      case cont_model::hill:
   
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = bkground; scaleMatrix(2,2)*= max_dose; 
        COV = scaleMatrix*COV*scaleMatrix.transpose(); 
        break; 
      case cont_model::exp_3:
       
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = 1/max_dose;  
        COV = scaleMatrix*COV*scaleMatrix.transpose(); 
        break; 
      case cont_model::exp_5:
      
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = 1/max_dose; 
        COV = scaleMatrix*COV*scaleMatrix.transpose(); 
        break; 
      case cont_model::power: 
      
        parms(0,0) *= bkground; parms(1,0) *= bkground*pow(1/max_dose,parms(2,0)); 
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = bkground*pow(1/max_dose,parms(2,0));
        scaleMatrix(1,2) = bkground*parms(1,0)*log(1/max_dose)*pow(1/max_dose,parms(2,0));  
        COV = scaleMatrix*COV*scaleMatrix.transpose(); 
        break; 
      case cont_model::polynomial:
      default:
        for (int i = 1; i <= degree; i++){
          scaleMatrix(i,i) *= pow(1/max_dose,i); 
        }
        COV = scaleMatrix*COV*scaleMatrix.transpose(); 
        break; 
    }
    return COV; 

  }
 

void rescale_mcmc(mcmcSamples *a, cont_model model,
                   double max_dose, bool is_logNormal,int degree){
   
   a->map_cov      = rescale_cov_matrix(a->map_cov, 
                                        a->map_estimate, model,
                                        max_dose, 1.0,
                                        is_logNormal, degree);
   a->map_estimate = rescale_parms(a->map_estimate,  model, 
                                   max_dose, 1.0,
                                   is_logNormal,  degree); 
   
   for (int i = 0; i < a->samples.cols(); i++){
     a->BMD(0,i) = a->BMD(0,i)*max_dose; 
     a->samples.col(i).array() = rescale_parms(a->samples.col(i),  model, 
                                                max_dose, 1.0,
                                                is_logNormal,  degree); 
   }
   return; 
}
 