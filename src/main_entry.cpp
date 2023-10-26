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


#include <iostream>
#include <string>
#include <vector>

#include "bmds_entry.h"

#include  <statmod.h>

#include <log_likelihoods.h>
#include <normal_likelihoods.h>
#include <normalModels.h>
#include <binomModels.h>
#include <IDPrior.h>

#include "bmd_calculate.h"

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
#include "list_r_conversion.h"
#include "owenst_asa076.h"

using namespace Rcpp;


#define MAX_PARMS 32 // Should never get close to this many!!!

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: owenst_fn
// purpose: takes input, which is assumed to be correct,
// and then computes the corresponding Owens T of those inputs
// output: function evaluation as a double
// [[Rcpp::export(".owenst_fn")]]
double owenst_fn(double x, double fx)
{
	return tfn(x, fx);
}

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_single_dichotomous
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis.
// output: BMD analysis with the model specified by NumericVector model
//
// [[Rcpp::export(".run_single_dichotomous")]]
List run_single_dichotomous(NumericVector model,
                            Eigen::MatrixXd data, Eigen::MatrixXd pr,
                            NumericVector options1, IntegerVector options2) 
{
  dichotomous_analysis Anal; 
  Anal.BMD_type =  (options1[0]==1)?eExtraRisk:eAddedRisk;
  Anal.BMR      =  options1[0]; 
  Anal.alpha    =  options1[1];
  Anal.parms    = pr.rows(); 
  Anal.model    = (dich_model)model[0]; 
  Anal.Y        = new double[data.rows()] ; 
  Anal.n_group  = new double[data.rows()] ; 
  Anal.doses    = new double[data.rows()] ; 
  Anal.prior    = new double[pr.cols()*pr.rows()];
  Anal.prior_cols = pr.cols(); 
  Anal.n          = data.rows(); 
  Anal.degree =   pr.rows()-1; 

  if (Anal.model == dich_model::d_multistage){
    Anal.degree = Anal.parms - 1; 
  }

  for (int i = 0; i < data.rows(); i++){
    Anal.Y[i] = data(i,1); 
    Anal.n_group[i] = data(i,2); 
  }
  
  for (int i = 0; i < data.rows(); i++){
    Anal.doses[i] = data(i,0); 
  }

  cp_prior(pr,Anal.prior);
  
  dichotomous_model_result res; 
  res.parms = new double[pr.rows()]; 
  res.cov   = new double[pr.rows()*pr.rows()]; 
  res.dist_numE = 200; 
  res.bmd_dist = new double[res.dist_numE*2]; 
  
  estimate_sm_laplace(&Anal, &res); 

  dichotomous_PGOF_data GOFdata; 
  GOFdata.n = Anal.n; GOFdata.Y = Anal.Y;
  GOFdata.model = Anal.model;
  GOFdata.model_df = res.model_df; 
  GOFdata.est_parms = res.parms; 
  GOFdata.doses = Anal.doses; GOFdata.n_group = Anal.n_group; 
  GOFdata.parms = Anal.parms; 
  dichotomous_PGOF_result GOFres; 
  
  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  GOFres.expected = new double[Anal.n]; 
  GOFres.residual = new double[Anal.n]; 

 
  compute_dichotomous_pearson_GOF(&GOFdata,&GOFres); 
  res.gof_p_value = GOFres.p_value;
  res.gof_chi_sqr_statistic = GOFres.test_statistic;

  
  Eigen::VectorXd resid(Anal.n); 
  Eigen::VectorXd expec(Anal.n);
  
  dichotomous_aod AOD;
  
  deviance_dichotomous(&Anal,&AOD);
  //cout << AOD.A1 << " : " << AOD.A2 << endl; 
  
  /*for (int i = 0; i < Anal.n; i++){
    cout << GOFres.expected[i] << " : " << GOFres.residual[i] << endl;   
  }*/
  
  delete[] GOFres.expected;
  delete[] GOFres.residual; 
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
 
  List rV = convert_dichotomous_fit_to_list(&res); 
  
  delete[] Anal.Y; 
  delete[] Anal.n_group; 
  delete[] Anal.doses; 
  delete[] Anal.prior; 
  delete[] res.parms;   
  delete[] res.cov;    
  delete[] res.bmd_dist;  
  return rV;
}

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_single_continuous
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis.
// output: BMD analysis with the model specified by NumericVector model
// [[Rcpp::export(".run_continuous_single")]]
List run_continuous_single(IntegerVector model, 
                           Eigen::MatrixXd Y, Eigen::MatrixXd X,
                           Eigen::MatrixXd prior, NumericVector options,
                           IntegerVector dist_type){
    
    bool   is_increasing = (bool)options[4]; 	double alpha = (double)options[3];
    double tail_p = (double)options[2]; 	double bmrf  = (double)options[1];
    int    riskType = (int)options[0];   
    unsigned int samples = (unsigned int) options[5];
    bool isFast = (bool) options[6]; 
    int transform =  options[8];
    
    
    //////////////////////////////////////////
    /// Set up the analysis
    ////////////////////////////////////////////////
    continuous_analysis anal; 
    anal.Y       =    new double[Y.rows()]; 
    anal.n       =    Y.rows(); 
    anal.n_group =    new double[Y.rows()]; 
    anal.sd      =    new double[Y.rows()]; 
    anal.doses   =    new double[Y.rows()]; 
    anal.model   =    (cont_model) model[0]; 
    anal.disttype     = dist_type[0]; 
    anal.isIncreasing = is_increasing; 
    anal.alpha        = 0.005; //alpha for analyses; 
    anal.BMD_type     = riskType; 
    anal.BMR          = bmrf; 
    anal.samples      = samples; 
    anal.tail_prob    = tail_p; 
    anal.suff_stat    = Y.cols()==3;
    anal.parms        = prior.rows();
    anal.prior_cols   = prior.cols(); 
    anal.degree       = 0;
    anal.transform_dose = transform; 
    anal.prior   = new double[prior.rows()*prior.cols()]; 
    cp_prior(prior,anal.prior);
    //
    // Check on the polynomial stuff
    //

    if (anal.model == cont_model::polynomial){
      // figure out the degree
      if (anal.disttype == distribution::normal ){
        anal.degree = anal.parms -2; 
      }else if (anal.disttype == distribution::normal_ncv){
        anal.degree = anal.parms -3; 
      }else{
        //throw an error! can'd do log-normal polynomial
        stop("Polynomial-Log-normal models are not allowed.\n Please choose normal or normal non-constant variance.");
      }
      
    }
   
   
    for (int i = 0; i < Y.rows(); i++){
      anal.Y[i] = Y(i,0); 
      anal.doses[i] = X(i,0); 
      if (Y.cols() == 3){ //sufficient statistics
         anal.n_group[i] = Y(i,1);
         anal.sd[i]      = Y(i,2); 
      }
    }
    ////////////////////////////////////
    //////////////////////////////////////////
    /// Set up the analysis
    ////////////////////////////////////////////////
    continuous_analysis anal2; 
    anal2.Y       =    new double[Y.rows()]; 
    anal2.n       =    Y.rows(); 
    anal2.n_group =    new double[Y.rows()]; 
    anal2.sd      =    new double[Y.rows()]; 
    anal2.doses   =    new double[Y.rows()]; 
    anal2.model   =    (cont_model) model[0]; 
    anal2.disttype     = dist_type[0]; 
    anal2.isIncreasing = is_increasing; 
    anal2.alpha        = 0.005; //alpha for analyses; 
    anal2.BMD_type     = riskType; 
    anal2.BMR          = bmrf; 
    anal2.samples      = samples; 
    anal2.tail_prob    = tail_p; 
    anal2.suff_stat    = Y.cols()==3;
    anal2.parms        = prior.rows();
    anal2.prior_cols   = prior.cols(); 
    anal2.degree       = 0;
    anal2.transform_dose = transform; 
    anal2.prior   = new double[prior.rows()*prior.cols()]; 
    cp_prior(prior,anal2.prior);
    //
    // Check on the polynomial stuff
    //

    if (anal2.model == cont_model::polynomial){
      // figure out the degree
      if (anal2.disttype == distribution::normal ){
        anal2.degree = anal2.parms -2; 
      }else if (anal2.disttype == distribution::normal_ncv){
        anal2.degree = anal2.parms -3; 
      }else{
        //throw an error! can'd do log-normal polynomial
        stop("Polynomial-Log-normal models are not allowed.\n Please choose normal or normal non-constant variance.");
      }
      
    }
   
   
    for (int i = 0; i < Y.rows(); i++){
      anal2.Y[i] = Y(i,0); 
      anal2.doses[i] = X(i,0); 
      if (Y.cols() == 3){ //sufficient statistics
         anal2.n_group[i] = Y(i,1);
         anal2.sd[i]      = Y(i,2); 
      }
    }
    
    
    ////////////////////////////////////

    ////////////////////////////////////
    continuous_model_result *result = new_continuous_model_result(anal.model,
                                                                  anal.parms,
                                                                  200); //have 200 equally spaced values
    ////////////////////////////////////
    continuous_deviance aod1; 
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        {
        estimate_sm_laplace(&anal,result,isFast);
        }

        #pragma omp section
        {
        
        if (anal.disttype == distribution::log_normal){
        
          estimate_log_normal_aod(&anal2,
                                  &aod1);
        
        }else{
          estimate_normal_aod(&anal2,
                                &aod1);
        }
        }
    }
    }
   continuous_expected_result exp_r; 
   exp_r.expected = new double[anal.n]; exp_r.n = anal.n; 
   exp_r.sd       = new double[anal.n];
   continuous_expectation(&anal, result,
                                &exp_r);
   
   NumericMatrix AOD(5,2);
   AOD(0,0) = aod1.A1; AOD(0,1) = aod1.N1; 
   AOD(1,0) = aod1.A2; AOD(1,1) = aod1.N2;
   AOD(2,0) = aod1.A3; AOD(2,1) = aod1.N3;
   AOD(3,0) = aod1.R;  AOD(3,1) = aod1.NR;
   AOD(4,0) = exp_r.like; AOD(4,1) = result->model_df;

    List rV = convert_continuous_fit_to_list(result); 

    rV.push_back(AOD,"Deviance");
    delete[] exp_r.expected; 
    delete[] exp_r.sd; 
    del_continuous_model_result(result); 
    del_continuous_analysis(anal);
    del_continuous_analysis(anal2);
    ////////////////////////////////////
    return rV; 
    
}
