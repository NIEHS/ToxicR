#pragma once

#ifndef lognormal_likelihoodsH
#define lognormal_likelihoodsH

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>

#endif

#include <gsl/gsl_randist.h>
#include <math.h>
#include "log_likelihoods.h"
#include "cmodeldefs.h"

/*class normalLL : public LL
 * This class defines a normal log-likelihood where the data 
 * are either normal sufficient statistics (i.e., sufficient_statistics = true)
 * or         the actual data              (i.e., sufficient_statistics = false)
 * 
 */
class lognormalLL: public LL {
public:
	lognormalLL(){
	}; 
	
	lognormalLL(Eigen::MatrixXd tY, Eigen::MatrixXd tX, bool SS) : LL(tY, tX) {
		sufficient_statistics = SS; // if it is a sufficient statistics model 
	};
	int    nParms() { return 2; }; 			// the model assumes one mean
										    // and constan variance
	bool isSuffStat(){
			return sufficient_statistics; 
	}; 
	
						    	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		double mean = theta(0,0); 
		Eigen::MatrixXd rV = Eigen::MatrixXd::Ones(Y.rows(),1)*mean; 
		return rV;  // mean is constant											
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		double var = theta(1,0); 
		Eigen::MatrixXd rV = Eigen::MatrixXd::Ones(Y.rows(),1)*exp(var); 
		return rV;  // variance is constant											
	};

	

	double negLogLikelihood(Eigen::MatrixXd theta) {
		// get the mean and variance for each dose group
		Eigen::MatrixXd mu = mean(theta); 
		Eigen::MatrixXd var = variance(theta); 
		Eigen::MatrixXd returnV = Y.col(0)*0; // make a log likelihood value
								                         // for each observation
								            
	     if (sufficient_statistics){
				// this assumes the correct geometric mean and geometric 
				// standard deviation are specified
				returnV = -Y.col(2).array()*Y.col(0).array() 
						    + Y.col(2).array()*log(1/sqrt(2.0*M_PI))
						    -(Y.col(2).array() / 2.0)*log(var.array()) 
						    - (1 / (2.0*var.array()))*((Y.col(2).array() - 1)*pow(Y.col(1).array(), 2) 
						    + Y.col(2).array()*pow(Y.col(0).array() - mu.array(), 2));
         }else{
				// log normal likelihood
				Eigen::MatrixXd sqerr = pow(log(Y.col(0).array())-mu.array(),2); 
				returnV = -log(Y.col(0).array())-(0.5)*log(2*M_PI*var.array())-(1/(2.0*var.array()))*sqerr.array(); 
		 }
		 double temp = -returnV.sum(); 
		 return temp; //isnan(temp)? -std::numeric_limits<double>::has_infinity:temp;
	}; 

	protected: 
		bool sufficient_statistics; 

};
#endif

