//File: DichMultistageBMD_NC.h
//Purpose: Creates a Dichotomous Multistage BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/5/2018
//Changes:
//
//
//
#include "log_likelihoods.h"
#include "binomModels.h"

#pragma once
#ifndef DichMultistageBMD_NCH
#define DichMultistageBMD_NCH

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif


#define MULTISTAGE_G(X) 1.0/(1.0 + exp(-X))
#define MULTISTAGE_POLY HILL_A(X) X
#define MULTISTAGE_B(X,D) X.bottomRows(D)

class dich_multistageNC : public binomialBMD {
private:
	int degree;
public:

	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		for (int i = 0; i < theta.rows(); i++)
			x[i] = theta(i, 0);
		return x;
	}

	// stub function not needed for multistage
	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		for (int i = 0; i < theta.rows(); i++)
			x[i] = theta(i, 0);
		return x;
	}

	// multistage dose response function
	dich_multistageNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int T) : binomialBMD(tY, tX) {
		degree = T;
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), degree);

		for (int i = 0; i < degree; i++) {
			newX.col(i) = pow(X.array(), i + 1);
		}
		X = newX;
		
	};
     
     virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
          Eigen::MatrixXd newX(D.rows(), degree);
          
          for (int i = 0; i < degree; i++) {
               newX.col(i) = pow(D.array(), i + 1);
          }
          
          return newX; 
          
     }
	int    nParms() { return degree + 1 ; }; // assumes n x 2 matrix degree + 1

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = MULTISTAGE_G(theta(0,0));
		Eigen::MatrixXd BETA = MULTISTAGE_B(theta,degree);
		Eigen::MatrixXd p(X.rows(), 1);
		p = -X * BETA;
		p = g+(1-g)*(1-p.array().exp());

		return   p;
	};

	Eigen::MatrixXd XgivenD(double d) {
		Eigen::MatrixXd rV(1, degree);
		for (int i = 0; i < degree; i++) {
			rV(0, i) = pow(d, i + 1);
		}
		// dose is the last parameter in X
		return rV;
	}

	bool fixedNC_BMDformula() { return false; } // for the multistage this is not true

												// Computes the EQUALITY BMD with a given BMR assuming no covariates NC option
	virtual double compute_BMD_ADDED_NC_EQUALITY(Eigen::MatrixXd theta, double *grad, double BMD, double BMR) {
		double g = MULTISTAGE_G(theta(0, 0));
		double tmp = exp(theta(0, 0));
		Eigen::MatrixXd BETA = MULTISTAGE_B(theta, degree);
		Eigen::MatrixXd doseZ = XgivenD(BMD);
		if (grad) {
			grad[0] = BMR * tmp / (BMR*tmp + BMR - 1);   // add extra
			for (int i = 0; i < degree; i++) {
				grad[i + 1] = doseZ(0, i);
			}
		}
		double rV = (doseZ*BETA)(0, 0);
		return rV + log(1 - BMR/(1-g));
	}

	virtual double compute_BMD_EXTRA_NC_EQUALITY(Eigen::MatrixXd theta, double *grad, double BMD, double BMR) {
		double g = MULTISTAGE_G(theta(0, 0));
		Eigen::MatrixXd BETA = MULTISTAGE_B(theta, degree);
		Eigen::MatrixXd doseZ = XgivenD(BMD);
		if (grad) {
			grad[0] = 0.0;
			for (int i = 0; i < degree; i++) {
				grad[i + 1] = doseZ(0, i);
			}
		}
		double rV = (doseZ*BETA)(0, 0);
		return rV + log(1-BMR);
	}


	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR , double isExtra) {
		double rV = 0.0;
		return rV;
	}

	// this is used for the starting value algorithm for the profile likelihood
	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		Eigen::MatrixXd temp = theta;
		double g = MULTISTAGE_G(theta(0, 0));
		double running_Val = 0.0;
		double A = BMR / (1.0 - g);
		// we modify the last values to get a BMD that works
		for (int i = 2; i < temp.rows(); i++) {
			running_Val += temp(i, 0)*pow(BMD, i);
		}
		temp(1, 0) = (-log(1.0-A)-running_Val)/BMD;
		return temp;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {

		Eigen::MatrixXd temp = theta;
		double g = MULTISTAGE_G(theta(0, 0));
		double running_Val = 0.0;
		double A = BMR;
		// we modify the last values to get a BMD that works
		for (int i = 2; i < temp.rows(); i++) {
			running_Val += temp(i, 0)*pow(BMD, i);
		}
		temp(1, 0) = (-log(1.0 - A) - running_Val) / BMD;
		return temp;
	}


	double compute_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double *grad) {

		double g = MULTISTAGE_G(theta(0, 0));
		Eigen::MatrixXd BETA = MULTISTAGE_B(theta, degree);

		double Z = 1.0;
		//grad is non null compute the gradient of the inequality
		// constraint
		if (grad) {
		}

		Z = Z / log(BMD);

		double rV = 0.0;
		if (geq) { // greater than or equal
			rV = inequality - Z + 1e-6;
			if (grad) {
			}
		}
		else {
			rV = Z - inequality + 1e-6;
			if (grad) {
			}

		}

		return rV;
	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		double g = MULTISTAGE_G(theta(0, 0));
		Eigen::MatrixXd BETA = MULTISTAGE_B(theta, degree);

		double Z = 1.0;
		Z = Z / log(BMD);
		double rV = 0.0;

		//grad is non null compute the gradient
		if (grad) {

		}

		if (geq) { // greater than or equal
			rV = inequality - Z + 1e-8;
			if (grad) {
				}
		}
		else {
			rV = Z - inequality + 1e-8;
			if (grad) {
			}

		}

		return rV;
	}

	int parameterRemoved() {
		return 3; // the Beta is the removed parameter for the Hill
	}


};

//Make sure the # defines are only for this file!
#endif
