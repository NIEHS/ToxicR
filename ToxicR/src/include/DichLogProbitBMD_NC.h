//File: DichLOGPROBITlBMD_NC.h
//Purpose: Creates a Dichotomous Hill BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/5/2018
//Changes:
//
//
//
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#define LOGPROBIT_G(X) 1.0/(1.0 + exp(-X))
#define LOGPROBIT_A(X) X
#define LOGPROBIT_B(X) X
#define LOGPROBIT_ADDED_Z(G,A,BMR) (gsl_cdf_gaussian_Pinv(BMR/(1-G), 1.0)-A)
#define LOGPROBIT_EXTRA_Z(G,A,BMR) (gsl_cdf_gaussian_Pinv(BMR, 1.0)-A)
#define LOGPROBIT_ADDED_RISK(G,A,B,BMR) exp(LOGPROBIT_ADDED_Z(G,A,BMR)/B)
#define LOGPROBIT_EXTRA_RISK(G,A,B,BMR) exp(LOGPROBIT_EXTRA_Z(G,A,BMR)/B)
#define LOGPROBIT_MEAN(G,A,B,D) ((D) <= (0.0)) ? (G) : (G+(1-G)*gsl_cdf_gaussian_P(A+B*log(D), 1))

#pragma once
#ifndef DichLogProbitBMD_NCH
#define DichLogProbitBMD_NCH

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif



#include "log_likelihoods.h"
#include "binomModels.h"

//void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd, void*)> math_func)
struct log_probit_inequality {
	double BMD;
	double BMR;
	bool geq;
	double inequality;
};

double logProbit_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);
double logProbit_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);

class dich_logProbitModelNC : public binomialBMD {

public:
	// hill dose response function
	dich_logProbitModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,  int T) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 3);
		Eigen::MatrixXd one(temp.rows(), 1);
		newX << one, one, temp;
		X = newX;
	};

	int    nParms() { return X.cols(); };

	virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
	     Eigen::MatrixXd returnV(D.rows(), 3);
	     Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1) ;
	     returnV << one , one , D.col(0); 
	     return returnV; 
	     
	}
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0)); double b = LOGPROBIT_B(theta(2, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = LOGPROBIT_MEAN(g, a, b, X(i, 2));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR, double isExtra) {
		double g = LOGPROBIT_G(theta(0, 0));
//		double a = LOGPROBIT_A(theta(1, 0));
//		double b = LOGPROBIT_B(theta(2, 0));


		double rV = (isExtra) ? -1 : BMR / (1 - g) -1;

		//if gradient is specified we calculate the gradient
		if (grad) {
			if (isExtra) {
				grad[0] = 0.0; grad[1] = 0.0;
			}
			else {
				grad[0] = -BMR * exp(theta(0, 0)) / pow(BMR + exp(theta(0, 0)), 2.0);
				grad[1] = 0.0;
			}
		}
		return rV;
	}

	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR) {
		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0));
		double b = LOGPROBIT_B(theta(2, 0));

		double BMD = LOGPROBIT_ADDED_RISK(g, a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0));
		double b = LOGPROBIT_B(theta(2, 0));

		double BMD = LOGPROBIT_EXTRA_RISK(g, a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }


	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0));
		double Z, temp;
		if (isExtra) {
			Z = LOGPROBIT_EXTRA_Z(g, a, BMR);
			temp = (lb*log(BMD) - Z) *1.2; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		else {
			Z = LOGPROBIT_ADDED_Z(g, a, BMR);
			temp = (lb*log(BMD) - Z) * 1.2; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0));
		double Z, temp;
		if (isExtra) {
			Z = LOGPROBIT_EXTRA_Z(g, a, BMR);
			temp = (Z - ub * log(BMD))*1.01; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		else {
			Z = LOGPROBIT_ADDED_Z(g, a, BMR);
			temp = (Z - ub * log(BMD))*1.01; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0) - temp;
		}
		return x;

	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0));


		double Z = LOGPROBIT_ADDED_Z(g, a, BMR);
		double BETA = Z / log(BMD);

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = theta(1, 0);
		newM(2, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {

		double g = LOGPROBIT_G(theta(0, 0));
		double a = LOGPROBIT_A(theta(1, 0));


		double Z = LOGPROBIT_EXTRA_Z(g, a, BMR);

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

		log_probit_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, logProbit_BMD_EXTRA_NC_INEQUALITY);
		}
		return logProbit_BMD_EXTRA_NC_INEQUALITY(theta, &M);

	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		log_probit_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, logProbit_BMD_ADDED_NC_INEQUALITY);
		}
		return logProbit_BMD_ADDED_NC_INEQUALITY(theta, &M);
	}

	int parameterRemoved() {
		return 2; // the Beta is the removed parameter for the LOGPROBIT
	}
};

//Make sure the # defines are only for this file!


#endif
