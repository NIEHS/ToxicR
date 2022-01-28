#pragma once

#ifndef normal_EXP_NCH
#define normal_EXP_NCH

#include "stdafx.h" // Precompiled header - does nothing if building R version
#include "normalModels.h"

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>

#define NORMAL_EXP2_UP 2
#define NORMAL_EXP2_DOWN 21
#define NORMAL_EXP3_UP 3
#define NORMAL_EXP3_DOWN 31
#define NORMAL_EXP4_UP 4
#define NORMAL_EXP4_DOWN 41
#define NORMAL_EXP5_UP 5
#define NORMAL_EXP5_DOWN 51

/* class normalEXPONENTIAL_BMD_NC : public LL
 *
 * This class defines a normal log-likelihood where 
 * the EXPONENTIAL model is the model that is fit to the data.
 * NC stands for no covariates, for future extentions that 
 * may include covariates.
 */

/*implements the basic model model Y = X*/
class normalEXPONENTIAL_BMD_NC : public normalLLModel {
	public: 

	normalEXPONENTIAL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS,
							            bool CV, int degree) : normalLLModel(tY, tX,SS,CV) {
			// if it is a sufficient statistics model 
			deg = degree; 
	};
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	}
	normalEXPONENTIAL_BMD_NC(){
	
	};
	
	virtual int parameter_to_remove(contbmd TYPE);

	normalEXPONENTIAL_BMD_NC(normalEXPONENTIAL_BMD_NC &M){
		 sufficient_statistics = M.sufficient_statistics; 
		 Y = M.Y; 
		 X = M.X; 
		 deg = M.deg; 
		 constant_variance = M.constant_variance; 
	};
	
	int  nParms() { 
				int add; 
				if (isConstVar()){
					add = 5; // Power model regression + constant variance
				} else{
					add = 6; // Hill regression + variance proportional to mean
				}; 		
			return add; 
	}
	
	virtual cont_model mean_type(){
	  switch(deg){
    	  case NORMAL_EXP3_DOWN:
    	  case NORMAL_EXP3_UP:
    	    return cont_model::exp_3;
    	    break; 
    	  case NORMAL_EXP5_DOWN:
    	  case NORMAL_EXP5_UP:
    	    return cont_model::exp_5; 
    	  default:
    	    return cont_model::generic;     
	  }
	  return cont_model::generic;
	}		
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta,Eigen::MatrixXd d){
				
		Eigen::MatrixXd rV; 
	  
		double sign = 1.0; 
		switch(deg){
			case NORMAL_EXP2_DOWN:
				sign = -1.0; 
			case NORMAL_EXP2_UP:
				rV = theta(0,0)*exp(sign*theta(1,0)*d.array()); 
				break; 
			case NORMAL_EXP3_DOWN:
				sign = -1.0;
			case NORMAL_EXP3_UP:
				rV = theta(0,0)*exp(sign*pow(theta(1, 0)*d.array(),theta(3,0)));
				break; 
			case NORMAL_EXP4_DOWN:
			case NORMAL_EXP4_UP:
				rV = theta(0,0)*(exp(theta(2,0))-(exp(theta(2,0))-1.0)*exp(-theta(1,0)*d.array()));
				break; 
			case NORMAL_EXP5_DOWN:
			case NORMAL_EXP5_UP:
			default:
				rV = theta(0,0)*(exp(theta(2,0))-(exp(theta(2,0))-1.0)*exp(-pow(theta(1,0)*d.array(),theta(3,0))));
			
		}
		return rV; 
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
		int deg; 
	// add BMD specific stuff here
}; 

#endif
