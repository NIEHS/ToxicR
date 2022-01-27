//File: log_likelihoods.h
//Purpose: Creates a set of log likelihoods that can be used
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 12/18/2017
//Changes:
//
// Date :05/06/18 Matt Wheeler
//       Added returnX function to the class LL
#pragma once

#ifndef log_likelihoodsH
#define log_likelihoodsH
#include "stdafx.h"
#include "gradient.h"

#include <cmath>
#include <nlopt.hpp>

#ifdef R_COMPILATION
    //necessary things to run in R
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else
    #include <Eigen/Dense>
    #include <gsl/gsl_randist.h>
#endif


class LL
{
public:
	LL(){};
	LL(Eigen::MatrixXd tY, Eigen::MatrixXd tX) : Y(tY), X(tX){};
	virtual int    nParms() { return 0; }; // return the number of parameters in the model
	int    nData() { return Y.rows(); }; // returns the number of datapoints in the model
	virtual double negLogLikelihood(Eigen::MatrixXd theta) { return 0.0; }; // no likelihood
     virtual Eigen::MatrixXd fixTheta(Eigen::MatrixXd theta){ return theta;}  // fix the parameters in the likelihood for some 
                                                                              // models; 

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return 0 * Y;  // no likelihood has no mean
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return 0 * Y;  // no likelihood has no variance
	};

	Eigen::MatrixXd      returnX(){
				return X;
	}
protected:
	Eigen::MatrixXd    Y; // Data - First column is observed, Second Column is N
	Eigen::MatrixXd		 X; // Covariates

};


class binomialLL: public LL {
public:
	binomialLL(Eigen::MatrixXd tY, Eigen::MatrixXd tX) : LL(tY, tX) {};
	int    nParms() { return Y.rows(); }; // the model is fully saturated

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return theta;  // mean is the sample mean at each dose
					   // group
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return theta.array()*(1-theta.array());  // no likelihood has no variance
	};


	double negLogLikelihood(Eigen::MatrixXd theta) {
		Eigen::MatrixXd p = mean(theta);
		Eigen::MatrixXd returnV = p;

		returnV = Y.col(0).array()*p.array().log() + (Y.col(1).array() - Y.col(0).array())*(1 - p.array()).log();

		// correct the likelihood 
		// "?log(0) is not there
		for (int i = 0; i < returnV.rows(); i++) {
          if (p(i, 0) < NEAR_ZERO)
            returnV(i, 0) = Y(i, 0)*log(NEAR_ZERO); // Note the last term doesn't exist because log(1) = 0
          else if ((1.0 - p(i, 0)) < NEAR_ZERO)
            returnV(i, 0) = (Y(i, 1) - Y(i, 0))*log(NEAR_ZERO); // Note the last term doesn't exist because log(1) = 0
        }

		return -returnV.sum();
	};

};
#endif
///////////////////////////////////////////////////////////////////////////////////
//
//
///////////////////////////////////////////////////////////////////////////////////
