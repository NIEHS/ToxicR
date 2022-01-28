#pragma once

#ifndef lognormal_HILL_NCH
#define lognormal_HILL_NCH

#include "log_likelihoods.h"
#include "lognormalModels.h"

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include "gsl/gsl_cdf.h"
#include <gsl/gsl_randist.h>

/*class lognormalHILL_BMD_NC : public LL
 * This class defines a normal log-likelihood where 
 * the HILL model is the model that is fit to the data.
 * NC stands for no covariates, for future extentions that 
 */
class lognormalHILL_BMD_NC : public lognormalLLModel {
	public: 

	// default constructor shouldn't ever be used
	lognormalHILL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, bool SS, int junk)
						        	: lognormalLLModel(tY, tX,SS) {
		  // if it is a sufficient statistics model 
	};
	
	lognormalHILL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, bool SS, bool CV, int junk)
	             : lognormalLLModel(tY, tX,SS) {
	  // if it is a sufficient statistics model 
	};
	
	virtual cont_model mean_type(){
	  return cont_model::hill; 
	}	
	
	lognormalHILL_BMD_NC(){
	
	};
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	}

	lognormalHILL_BMD_NC(lognormalHILL_BMD_NC &M){
		 sufficient_statistics = M.sufficient_statistics; 
		 Y = M.Y; 
		 X = M.X; 
	};
	
	int    nParms() { 
			return 5; // Hill regression + constant variance
	}
	
	virtual int parameter_to_remove(contbmd TYPE);

	virtual Eigen::MatrixXd mmean(Eigen::MatrixXd theta,Eigen::MatrixXd d){
	  double gamma = theta(0,0); 
	  double nu = theta(1,0); 
	  double k  = theta(2,0); 
	  double n_exp = theta(3,0); 
	  double var_t  = exp(theta(4,0)); 
	  
	  Eigen::MatrixXd rV = gamma + nu*pow(d.array(),n_exp)/(pow(k,n_exp)+pow(d.array(),n_exp));		
	  
	  return exp(log(rV.array())+0.5*var_t); 
	}
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta,Eigen::MatrixXd d){
		double gamma = theta(0,0); 
		double nu = theta(1,0); 
		double k  = theta(2,0); 
		double n_exp = theta(3,0); 
		
		Eigen::MatrixXd rV = gamma + nu*pow(d.array(),n_exp)/(pow(k,n_exp)+pow(d.array(),n_exp));		
		
		return log(rV.array()); 
	}

	// return true if it is a increasing function

	//BASIC BMD Computation absolute to hybrid
	double bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing); 
	double bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing); 
	double bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing); 
	double bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing); 
	double bmd_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing); 
    double bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,double BPROB);
    
    //BASIC BMD Computation absolute to hybrid
	double bmd_absolute_bound(Eigen::MatrixXd theta,	 double BMD,double BMRF, bool isIncreasing); 
	double bmd_stdev_bound(Eigen::MatrixXd theta, 	     double BMD,double BMRF, bool isIncreasing); 
	double bmd_reldev_bound(Eigen::MatrixXd theta,       double BMD,double BMRF, bool isIncreasing); 
	double bmd_point_bound(Eigen::MatrixXd theta,        double BMD,double BMRF, bool isIncreasing); 
	double bmd_extra_bound(Eigen::MatrixXd theta,        double BMD,double BMRF, bool isIncreasing); 
    double bmd_hybrid_extra_bound(Eigen::MatrixXd theta, double BMD,double BMRF, bool isIncreasing,
								  double TAIL_PROB);

	///
	/*Eigen::MatrixXd bmd_start_extra_hybrid(Eigen::MatrixXd theta, double BMD,
											double BMRF, bool isIncreasing,
											double TAIL_PROB);*/
	//
	virtual double bmd_start_absolute(unsigned n,
									  const double *b,
									  double *grad,
									  void   *data); 
	virtual std::vector<double>  bmd_start_absolute_clean(std::vector<double> x,
															double BMRF, double BMD, bool isIncreasing);
	
	//
	/*Eigen::MatrixXd bmd_start_stdev(Eigen::MatrixXd theta, double BMD,
									double BMRF, bool isIncreasing,
									double TAIL_PROB);
	//*/
	virtual double bmd_start_reldev(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_reldev_clean(std::vector<double> x,
										  double BMRF, 
										  double BMD, 
										  bool isIncreasing);

	/////////////////////////////////////////////////////////////////////////
	// standard deviation starting value code
	//
	/////////////////////////////////////////////////////////////////////////
	virtual double bmd_start_stddev(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_stddev_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing);

	/////////////////////////////////////////////////////////////////////////
	// standard deviation starting value code
	/////////////////////////////////////////////////////////////////////////
	
	virtual double bmd_start_extra(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_extra_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing);

	/////////////////////////////////////////////////////////
	virtual double bmd_start_point(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_point_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing);

	/////////////////////////////////////////////////////////
	
	virtual double  bmd_start_hybrid_extra(unsigned n,
													  const double *b,
													  double *grad,
													  void   *data); 
													  
	std::vector<double> bmd_start_hybrid_extra_clean(std::vector<double> x,
						  							 double BMRF, double BMD,
													 bool isIncreasing, double tail_prob);												  

	/////////////////////////////////////////////////////////
	/*
	Eigen::MatrixXd bmd_start_point(Eigen::MatrixXd theta, double BMD,
									double BMRF, bool isIncreasing,
									double TAIL_PROB);
	//
	Eigen::MatrixXd bmd_start_extra(Eigen::MatrixXd theta, double BMD,
									double BMRF, bool isIncreasing,
									double TAIL_PROB);*/

	private:
	// add BMD specific stuff here
}; 

#endif

