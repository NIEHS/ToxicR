#include "log_likelihoods.h"


#pragma once
#ifndef normal_likelihoodsH
#define normal_likelihoodsH


#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif


#include <gsl/gsl_randist.h>
#include "cmodeldefs.h"

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

/*class normalLL : public LL
 * This class defines a normal log-likelihood where the data 
 * are either normal sufficient statistics (i.e., sufficient_statistics = true)
 * or         the actual data              (i.e., sufficient_statistics = false)
 * 
 */

class normalLL: public LL {
public:
	normalLL(){}; 

	normalLL(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS) : LL(tY, tX) {
		sufficient_statistics = SS; // if it is a sufficient statistics model 
 
	};
	
	
	
	int    nParms() { return 2; }; 		// the model assumes one mean
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
		Eigen::MatrixXd rV = Eigen::MatrixXd::Ones(Y.rows(),1)*var; 
		return rV;  // variance is constant											
	};
  
	double negLogLikelihood(Eigen::MatrixXd theta);
	protected: 
		bool sufficient_statistics; 


};
#endif
