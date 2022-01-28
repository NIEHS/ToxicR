//File: DichWEIBULLlBMD_NC.h
//Purpose: Creates a Dichotomous Hill BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/5/2018
//Changes:
//
//
//
#define WEIBULL_G(X) 1.0/(1.0 + exp(-X))
#define WEIBULL_A(X) X
#define WEIBULL_B(X) X
#define WEIBULL_EXTRA_Z(G,A,BMR) pow(-log(1.0-BMR),1.0/A)
#define WEIBULL_ADDED_Z(G,A,BMR) pow(-log(1.0-BMR/(1-G)),1.0/A)
#define WEIBULL_ADDED_RISK(G,A,B,BMR) WEIBULL_ADDED_Z(G,A,BMR)/pow(B,1.0/A)
#define WEIBULL_EXTRA_RISK(G,A,B,BMR) WEIBULL_EXTRA_Z(G,A,BMR)/pow(B,1.0/A)
#define WEIBULL_MEAN(G,A,B,D) ((D) <= (0.0)) ? (G) : (G+(1-G)*(1-exp(-B*pow(D,A))))

#pragma once
#ifndef DichWeibullBMD_NCH
#define DichWeibullBMD_NCH

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
struct log_weibull_inequality {
	double BMD;
	double BMR;
	bool geq;
	double inequality;
};

double WEIBULL_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);

double WEIBULL_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);

class dich_weibullModelNC : public binomialBMD {

public:
	// hill dose response function
	dich_weibullModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,int Deg) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 3);
		Eigen::MatrixXd one(temp.rows(), 1);
		newX << one, one, temp;
		X = newX;
	};
	int    nParms() { return X.cols(); };

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
	     Eigen::MatrixXd returnV(D.rows(), 3);
	     Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1) ;
	     returnV << one , one , D.col(0); 
	     return returnV; 
	     
	}
	
	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0)); double b = WEIBULL_B(theta(2, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = WEIBULL_MEAN(g, a, b, X(i, 2));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR, double isExtra) {
		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0));
		double b = WEIBULL_B(theta(2, 0));


		double rV = (isExtra) ? -1.0 : BMR / (1 - g) -1;

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
		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0));
		double b = WEIBULL_B(theta(2, 0));

		double BMD = WEIBULL_ADDED_RISK(g, a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0));
		double b = WEIBULL_B(theta(2, 0));

		double BMD = WEIBULL_EXTRA_RISK(g, a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }


	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0));
		double Z, temp;
		if (isExtra) {
			Z = pow(WEIBULL_EXTRA_Z(g, a, BMR),a);
			temp = Z/pow(BMD,a); // temp is the 'new beta with constraint satisfied
			if (temp < lb){
					temp = Z/pow(BMD,a) + fabs( Z/pow(BMD,a) - lb)*1.2;
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}else{
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}
		}
		else {
			Z = pow(WEIBULL_ADDED_Z(g, a, BMR),a);
			// lco: Added the line below to resolve using "temp" before it is initialized.
			// It was not in Matt's latest code as of 7/2/2018.
			temp = Z / pow(BMD, a); // temp is the 'new beta with constraint satisfied
			if (temp < lb){
					temp = Z/pow(BMD,a) + fabs( Z/pow(BMD,a) - lb)*1.2;
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}else{
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}
		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		double g = WEIBULL_G(theta(0, 0));
		double a = 1;//WEIBULL_A(theta(1, 0));
		double Z, temp;
		if (isExtra) {
			Z = pow(WEIBULL_EXTRA_Z(g, a, BMR),a);
			temp = Z/pow(BMD,a); // temp is the 'new beta with constraint satisfied
			if (temp > ub){
					temp = Z/pow(BMD,a) - fabs( Z/pow(BMD,a) - ub)*1.2;
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}else{
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}
		}
		else {
			Z = pow(WEIBULL_ADDED_Z(g, a, BMR),a);
			temp = Z/pow(BMD,a); // temp is the 'new beta with constraint satisfied
			if (temp > ub){
					temp = Z/pow(BMD,a) - fabs( Z/pow(BMD,a) - ub)*1.2;
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}else{
					x[0] = theta(0, 0);
					x[1] = -(log(temp) -log(Z))/log(BMD);
			}
		}
		return x;

	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0));


		double Z = WEIBULL_ADDED_Z(g, a, BMR);
		double BETA = pow(Z,a) / pow(BMD,a);

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = theta(1, 0);
		newM(2, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {

		double g = WEIBULL_G(theta(0, 0));
		double a = WEIBULL_A(theta(1, 0));


		double Z = WEIBULL_EXTRA_Z(g, a, BMR);

		double BETA =  pow(Z, a) / pow(BMD, a);

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

		log_weibull_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, WEIBULL_BMD_EXTRA_NC_INEQUALITY);
		}

		return WEIBULL_BMD_EXTRA_NC_INEQUALITY(theta, &M);

	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		log_weibull_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, WEIBULL_BMD_ADDED_NC_INEQUALITY);
		}
		return WEIBULL_BMD_ADDED_NC_INEQUALITY(theta, &M);
	}

	int parameterRemoved() {
		return 2; // the Beta is the removed parameter for the WEIBULL
	}
};

//Make sure the # defines are only for this file!

#endif
