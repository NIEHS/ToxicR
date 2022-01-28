
#ifndef CONTINUOUS_CLEAN_AUX_H
#define CONTINUOUS_CLEAN_AUX_H
#include <vector>

#ifdef R_COMPILATION  
  #include <RcppEigen.h>
  #include <RcppGSL.h>
#else 
  #include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>
#include "bmdStruct.h"
#include "mcmc_struct.h"

std::vector<double> unique_list(Eigen::MatrixXd X); 
Eigen::MatrixXd cleanSuffStat(Eigen::MatrixXd Y, Eigen::MatrixXd X, bool is_logNormal,bool use_divisor = true); 

double get_divisor(Eigen::MatrixXd Y, Eigen::MatrixXd X);
Eigen::MatrixXd createSuffStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                               bool is_logNormal);

Eigen::MatrixXd rescale_parms(Eigen::MatrixXd parms, cont_model model,
                              double max_dose, double bkground,bool is_logNormal,
                              int degree = 2); 
                              
Eigen::MatrixXd rescale_cov_matrix(Eigen::MatrixXd COV, 
                									Eigen::MatrixXd parms, cont_model model,
                									double max_dose, double bkground,
                									bool is_logNormal, int degree = 2);

void rescale_mcmc(mcmcSamples *a, cont_model model,
                  double max_dose, bool is_logNormal,int degree = 2); 


#endif
