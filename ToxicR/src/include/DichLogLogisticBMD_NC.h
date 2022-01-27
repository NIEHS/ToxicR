//File: DichLogLogisticlBMD_NC.h
//Purpose: Creates a Dichotomous Hill BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/5/2018
//Changes:
//
//
//
#pragma once
#ifndef DichLogLogisticBMD_NCH
#define DichLogLogisticBMD_NCH

#include <iostream>

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif

#include "log_likelihoods.h"
#include "binomModels.h"

#define LOGLOGISTIC_G(X) 1.0/(1.0 + exp(-X))
#define LOGLOGISTIC_A(X) X
#define LOGLOGISTIC_B(X) X
#define LOGLOGISTIC_ADDED_Z(G,A,BMR) (log(BMR/(1-G-BMR))-A)
#define LOGLOGISTIC_EXTRA_Z(G,A,BMR)  (log(BMR/(1-BMR))-A)
#define LOGLOGISTIC_ADDED_RISK(G,A,B,BMR) exp(LOGLOGISTIC_ADDED_Z(G,A,BMR)/B)
#define LOGLOGISTIC_EXTRA_RISK(G,A,B,BMR) exp(LOGLOGISTIC_EXTRA_Z(G,A,BMR)/B)
#define LOGLOGISTIC_MEAN(G,A,B,D) ((D) <= (0.0)) ? (G) : (G+(1-G)/(1+exp(-A-B*log(D))))


class dich_loglogisticModelNC : public binomialBMD {

public:
	// hill dose response function
	dich_loglogisticModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int T) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 3);
		Eigen::MatrixXd one(temp.rows(), 1);
		newX << one, one, temp;
		X = newX;
	};
	int    nParms() { return X.cols();};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
	     Eigen::MatrixXd returnV(D.rows(), 3);
	     Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1) ;
	     returnV << one , one , D.col(0); 
	     return returnV; 
	     
	}
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1,0)); double b = LOGLOGISTIC_B(theta(2, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = LOGLOGISTIC_MEAN(g, a, b, X(i, 2));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR , double isExtra) {
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));
		double b = LOGLOGISTIC_B(theta(2, 0));


		double rV = (isExtra) ? -1.0 : -BMR / (1 - g - BMR);

		//if gradient is specified we calculate the gradient
		if (grad) {
			if (isExtra) {
				grad[0] = 0.0; grad[1] = 0.0;
			}
			else {
				grad[0] = -BMR*exp(theta(0,0))/pow(BMR+exp(theta(0, 0))-1.0,2.0);
				grad[1] = 0.0;
			}
		}
		return rV;
	}

	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR) {
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));
		double b = LOGLOGISTIC_B(theta(2, 0));

		double BMD = LOGLOGISTIC_ADDED_RISK(g, a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));
		double b = LOGLOGISTIC_B(theta(2, 0));

		double BMD = LOGLOGISTIC_EXTRA_RISK(g, a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }


	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR,bool isExtra,double lb) {
		std::vector<double> x(theta.rows());
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));
		double Z, temp;
		if (isExtra) {
			Z = LOGLOGISTIC_EXTRA_Z(g, a, BMR);
			temp = (lb*log(BMD) - Z)*2; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}else {
			Z = LOGLOGISTIC_ADDED_Z(g, a, BMR);
			temp = (lb*log(BMD) - Z)*2; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR,bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));
		double Z, temp;
		if (isExtra) {
			Z = LOGLOGISTIC_EXTRA_Z(g, a, BMR);
			temp = (Z-ub*log(BMD))*1.01; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		else {
			Z = LOGLOGISTIC_ADDED_Z(g, a, BMR);
			temp = (Z - ub*log(BMD))*1.01; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		return x;

	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));


		double Z = LOGLOGISTIC_ADDED_Z(g, a, BMR);
		double BETA = Z / log(BMD);

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = theta(1, 0);
		newM(2, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {

		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));


		double Z = LOGLOGISTIC_EXTRA_Z(g, a, BMR);

		double BETA = Z / log(BMD);

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = theta(1, 0);
		newM(2, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd XgivenD(double d) {
		Eigen::MatrixXd rV(1, 3); rV << 1.0, 1.0, d; // by default the
													 // dose is the last parameter in X
		return rV;
	}

	double compute_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double *grad) {

		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));

		double Z = LOGLOGISTIC_EXTRA_Z(g, a, BMR); //note BMD is a placeholder

		//grad is non null compute the gradient of the inequality
		// constraint
		if (grad) {
			grad[0] = 0; //
			grad[1] = -1.0; // alpha
		}

		Z = Z / log(BMD);

		double rV = 0.0;
		if (geq) { // greater than or equal
			rV = inequality - Z ;
			if (grad) {
				grad[0] = grad[0] * (-1.0 / log(BMD));
				grad[1] = grad[1] * (-1.0 / log(BMD));
			}
		}
		else {
			rV = Z - inequality;
			if (grad) {
				grad[0] = grad[0] * (1.0 / log(BMD));
				grad[1] = grad[1] * (1.0 / log(BMD));
			}

		}

		return rV;
	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		double g = LOGLOGISTIC_G(theta(0, 0));
		double a = LOGLOGISTIC_A(theta(1, 0));


		double Z = LOGLOGISTIC_ADDED_Z(g, a, BMR);
		Z = Z / log(BMD);
		double rV = 0.0;

		//grad is non null compute the gradient
		if (grad) {
			grad[0] = - exp(theta(0, 0)) /((1.0 + exp(theta(0, 0)))*(exp(theta(0, 0)) + 1.0));
			grad[1] = - 1.0; // alpha
		}

		if (geq) { // greater than or equal
			rV = inequality - Z;
			if (grad) {
				grad[0] *= (-1.0 / log(BMD));
				grad[1] *= (-1.0 / log(BMD));
			}
		}
		else {
			rV = Z - inequality ;
			if (grad) {
				grad[0] *= (1.0 / log(BMD));
				grad[1] *= (1.0 / log(BMD));
			}

		}

		return rV;
	}

	int parameterRemoved() {
		return 2; // the Beta is the removed parameter for the LogLogistic
	}
};

//Make sure the # defines are only for this file!
#undef LOGLOGISTIC_G
#undef LOGLOGISTIC_A
#undef LOGLOGISTIC_B
#undef LOGLOGISTIC_EXTRA_Z
#undef LOGLOGISTIC_ADDED_Z
#undef LOGLOGISTIC_ADDED_RISK
#undef LOGLOGISTIC_EXTRA_RISK
#endif
