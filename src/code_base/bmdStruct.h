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

//necessary things to run in R  
#ifndef BMD_ANALYSIS_h
#define BMD_ANALYSIS_h

#include "cmodeldefs.h"

// Dichotomous Structures
//
// dichotomous_analysis:   
//   Purpose - Contains all of the information for a dichotomous analysis.
//   It is used do describe a single model analysis, in which all of the information
//   is used, or a MA analysis, in which all the information save prior, degree, parms
//   and prior_cols are used. 
struct dichotomous_analysis{
  int model; // Model Type as listed in dich_model
  int n;     // total number of observations obs/n 
  double *Y; // observed +
  double *doses; // 
  double *n_group; //size of the group
  double *prior; // a column order matrix parms X prior_cols 
  int BMD_type; // 1 = extra ; added otherwise
  double BMR; 
  double alpha; // alpha of the analysis
  int degree;  // degree of polynomial used only  multistage
  int samples; // number of MCMC samples. 
  int burnin;  // size of burin
  int parms;   // number of parameters in the model
  int prior_cols; // colunns in the prior
};

struct dichotomous_PGOF_data{
     int     n;     // total number of observations obs/n 
     double *Y; // observed +
     double *doses; // 
     double *n_group; //size of the group
     double  model_df;
     int     model; 
     int     parms; 
     double  *est_parms; 
};

struct dichotomous_PGOF_result{
     int     n;        // total number of observations obs/n 
     double *expected; // 
     double *residual; //size of the group
     double  test_statistic; 
     double  p_value; 
     double  df;  
};

struct dichotomous_aod{
  double A1; 
  int     N1;        // total number of observations obs/n 
  double A2; 
  int N2;   
};

struct continuous_expected_result{
  int     n;        // total number of observations obs/n 
  double *expected; // 
  double *sd; 
  double like; 
};

//
// dichotomous_model_result: 
// Purpose: Data structure that is populated with all of the necessary
// information for a single model fit. 
//
struct dichotomous_model_result{
  int      model;               // dichotomous model specification
  int      nparms; 		          //number of parameters in the model
  double  *parms;               // Parameter Estimate 
  double  *cov;                 // Covariance Estimate
  double   max;                 // Value of the Likelihood/Posterior at the maximum
  int      dist_numE;           // number of entries in rows for the bmd_dist
  double      model_df;         // Used model degrees of freedom
  double      total_df;         // Total degrees of freedom
  double  *bmd_dist;            // bmd distribution (dist_numE x 2) matrix
  double  bmd;                  // the central estimate of the BMD
  double gof_p_value;           // P-value from Chi Square goodness of fit
  double gof_chi_sqr_statistic; // Chi Square Statistic for goodness of fit
};

//
// dichotomousMA_analysis
// Puprose: Fill out all of the information for a dichotomous
// model average.  
//
//
struct dichotomousMA_analysis{
  int    nmodels;      //number of models for the model average
  double **priors;     // List of pointers to prior arrays
                       // priors[i] is the prior array for the ith model ect
  int    *nparms;      //parameters in each model
  int    *actual_parms;//actual number of parameters in the model
  int    *prior_cols;  // columns in the prior if there are 'more' in the future
                       // presently there are only 5
  int    *models;      // list of models this is defined by dich_model. 
  double *modelPriors; // prior probability on the model
};

// 
//
//
struct dichotomousMA_result{
  int                       nmodels; //number of models for each 
  struct dichotomous_model_result **models; // Individual model fits for each
                                     // model average  
  int                     dist_numE; // number of entries in rows for the bmd_dist
  double                *post_probs; // posterior probabilities
  double                  *bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
};

// Continuous Structures
struct continuous_analysis{
  enum cont_model model; 
  int n; 
  bool suff_stat; //true if the data are in sufficient statistics format
  double *Y; // observed data means or actual data
  double *doses; // 
  double *sd; // SD of the group if suff_stat = true, null otherwise. 
  double *n_group; // N for each group if suff_stat = true, null otherwise
  double *prior; // a column order matrix px5 where p is 
                // the number of parametersd
  int BMD_type; // type of BMD
  
  bool isIncreasing; // if the BMD is defined increasing or decreasing
  double BMR; // Benchmark response related to the BMD type
  double tail_prob; // tail probability 
  int    disttype;  // Distribution type defined in the enum distribution
  double alpha;     // specified alpha
  int samples; // number of MCMC samples.
  int degree; // if polynomial it is the degree 
  int burnin;  // burn in 
  int parms; // number of parameters 
  int prior_cols; 
  int transform_dose; // Use the arc-sin-hyperbolic inverse to transform dose. 
};

struct continuousMA_analysis{
  int    nmodels;         //number of models for each 
  double **priors;        // pointer to pointer arrays for the prior
                          // each prior will have nparms[i] x prior_cols[i]
                          // it is also a columnwise matrix
  int    *nparms;         //parameter in each model
  int    *actual_parms;//actual number of parameters in the model
  int    *prior_cols;  // columns in the prior if there are 'more' in the future
  int    *models;      // given model
  int    *disttype;    // given distribution type
  double *modelPriors; // prior probability on the model
};

struct continuous_model_result{
  int      model;           // continuous model specification
  int      dist;            // distribution_type 
  int      nparms; 		      //number of parameters in the model
  double  *parms;           // Parameter Estimate 
  double  *cov;             // Covariance Estimate
  double   max;             // Value of the Likelihood/Posterior at the maximum
  int      dist_numE;       // number of entries in rows for the bmd_dist
  double    model_df;        // Used model degrees of freedom
  double    total_df;        // Total degrees of freedom
  double    bmd;             // The bmd at the maximum
  double   *bmd_dist;        // bmd distribution (dist_numE x 2) matrix
};

struct continuousMA_result{
  int                      nmodels; //number of models for each 
  struct continuous_model_result **models; //priors
  int                    dist_numE; // number of entries in rows for the bmd_dist
  double               *post_probs; // posterior probabilities
  double                 *bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
};

// mcmc structures
struct bmd_analysis_MCMC{
  int model; // model used in the analysis
  unsigned int burnin; // burnin samples
  unsigned int samples; // total samples including burnin
  unsigned int nparms;  // parameters in the model
  double * BMDS;        // array of samples of BMDS length (samples)
  double * parms;       // array of parameters length (samples X parms)
                        
};

struct ma_MCMCfits{
  unsigned int nfits; 
  struct bmd_analysis_MCMC **analyses;   
}; 

struct continuous_deviance{
  double A1; 
  int N1; 
  double A2; 
  int N2; 
  double A3; 
  int N3; 
  double R; 
  int NR; 
};


// odds and ends
struct bmd_analysis_MCMC * new_mcmc_analysis(int model,
                                       int parms, 
                                       unsigned int samples);
void del_mcmc_analysis(struct bmd_analysis_MCMC *an); 

#if defined(__cplusplus) || defined(R_COMPILATION)
void cp_prior(Eigen::MatrixXd temp, double *priors);
#endif

void del_continuousMA_analysis(struct continuousMA_analysis CMA);
void del_continuous_analysis(struct continuous_analysis a); 

struct dichotomousMA_result * new_dichotomousMA_result(int nmodels,
                                                int dist_numE);

struct continuous_model_result * new_continuous_model_result(int model,
													  unsigned int n_parm,
                            unsigned int n_elm);
                                                      
struct dichotomous_model_result * new_dichotomous_model_result(int model,
                                                        int parms,
                                                        int dist_numE);

void del_continuous_model_result(struct continuous_model_result * cm); 
void delete_dichotomousMA_result(struct dichotomousMA_result *res); 


#endif
