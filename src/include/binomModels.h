//File: log_likelihoods.h
//Purpose: Creates a set of log likelihoods that can be used
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 12/18/2017
//Changes:
//
//
//
#pragma once
#include "log_likelihoods.h"
#include "stdafx.h"

#ifndef binomialBMDH
#define binomialBMDH

#ifdef R_COMPILATION
    //necessary things to run in R
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else
    #include <Eigen/Dense>
#endif


class binomialBMD : public binomialLL {
public:
	binomialBMD(Eigen::MatrixXd tY, Eigen::MatrixXd tX) : binomialLL(tY, tX) {};

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	
	
	virtual bool fixedNC_BMDformula() { return false; }

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd X,double BMR) {
		return 0.0;
	}
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return variance(theta, X);
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
		Eigen::MatrixXd p = mean(theta); 
		return  p.array()*(1 - p.array());
	};
     virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
          Eigen::MatrixXd returnV(D.rows(), 2);
          Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1);
          returnV << one , D.col(0); 
          return returnV; 
          
     }
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
		return   1 / (1 + exp(-(X*theta).array()));
	};


	//computes an X given D
	virtual Eigen::MatrixXd XgivenD(double d) {
		Eigen::MatrixXd rV(1, 2); rV << 1.0, d; // by default the
												// dose is the last parameter in X
		return rV;
	} //TODO Make one with covariates!

	// Computes the BMD with a given BMR assuming no covariates NC option
	virtual double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR) {
		double tempD = 1.0; //start the BMD at MAX DOSE
		Eigen::MatrixXd doseZ = XgivenD(0.0);
		Eigen::MatrixXd dose = XgivenD(tempD);

		while ((mean(theta, dose)(0, 0) - mean(theta, doseZ)(0, 0)) < BMR) {
			if (tempD > 128.0) return tempD; //BMD is 128 times MAX DOSE!
			tempD *= 2.0;
			dose = XgivenD(tempD);
		}

		double max = tempD; double min = 0.0; double mpoint = (max + min) / 2.0;
		dose = XgivenD(mpoint);
		double test = mean(theta, dose)(0, 0) - mean(theta, doseZ)(0, 0) - BMR;
		while (fabs(test) > NEAR_ZERO) {
			if (test > 0.0)
				max = mpoint;
			else
				min = mpoint;

			mpoint = (max + min) / 2.0;
			dose = XgivenD(mpoint);
			test = mean(theta, dose)(0, 0) - mean(theta, doseZ)(0, 0) - BMR;
		}

		return mpoint;
	}

	virtual double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double tempD = 1.0; //start the BMD at MAX DOSE
		Eigen::MatrixXd doseZ = XgivenD(0.0);
		Eigen::MatrixXd dose = XgivenD(tempD);
		double a = mean(theta, dose)(0, 0);
		double b = mean(theta, doseZ)(0, 0);
		while (((a - b) / (1 - b)) < BMR) {
			if (tempD > 128.0) return tempD; //BMD is 128 times MAX DOSE!
			tempD *= 2.0;
			dose = XgivenD(tempD);
			a = mean(theta, dose)(0, 0);

		}

		double max = tempD; double min = 0.0; double mpoint = (max + min) / 2.0;
		dose = XgivenD(mpoint);
		a = mean(theta, dose)(0, 0);
		double test = (a - b) / (1 - b) - BMR;
		while (fabs(test) > NEAR_ZERO) {
			if (test > 0.0)
				max = mpoint;
			else
				min = mpoint;

			mpoint = (max + min) / 2.0;
			dose = XgivenD(mpoint);
			a = mean(theta, dose)(0, 0);
			test = (a - b) / (1 - b) - BMR;
		}

		return mpoint;
	}

	// Computes the EQUALITY BMD with a given BMR assuming no covariates NC option
	virtual double compute_BMD_ADDED_NC_EQUALITY(Eigen::MatrixXd theta,double *grad, double BMD, double BMR) {
		double tempD = 1.0; //start the BMD at MAX DOSE
		Eigen::MatrixXd doseZ = XgivenD(0.0);
		Eigen::MatrixXd dose = XgivenD(tempD);
		return ((mean(theta, dose)(0, 0) - mean(theta, doseZ)(0, 0)) - BMR);
	}

	virtual double compute_BMD_EXTRA_NC_EQUALITY(Eigen::MatrixXd theta, double* grad,double BMD, double BMR) {
		double tempD = 1.0; //start the BMD at MAX DOSE
		Eigen::MatrixXd doseZ = XgivenD(0.0);
		Eigen::MatrixXd dose = XgivenD(tempD);
		double a = mean(theta, dose)(0, 0);
		double b = mean(theta, doseZ)(0, 0);
		return (((a - b) / (1 - b)) - BMR);
	}

	};

#endif

