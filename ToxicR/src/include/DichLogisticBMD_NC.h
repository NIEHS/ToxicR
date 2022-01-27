//File: DichLOGISTIClBMD_NC.h
//Purpose: Creates a Dichotomous QuantalLinear BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/25/2018
//Changes:
//
//
//

#undef ADDED_Q
#undef EXTRA_Q
#define LOGISTIC_A(X) X
#define LOGISTIC_B(X) X
#define ADDED_Q(A,BMR) (BMR*(1.0+exp(-A))/exp(-A))
#define EXTRA_Q(A,BMR) (BMR)
#define LOGISTIC_ADDED_Z(A,BMR) -log((1-ADDED_Q(A,BMR))/(1+ADDED_Q(A,BMR)*exp(-A)))
#define LOGISTIC_EXTRA_Z(A,BMR) -log((1-EXTRA_Q(A,BMR))/(1+EXTRA_Q(A,BMR)*exp(-A)))
#define LOGISTIC_ADDED_RISK(A,B,BMR) LOGISTIC_ADDED_Z(A,BMR)/B
#define LOGISTIC_EXTRA_RISK(A,B,BMR) LOGISTIC_EXTRA_Z(A,BMR)/B
#define LOGISTIC_MEAN(A,B,D) 1/(1+exp(-A-B*D))

#pragma once
#ifndef DichLogisticBMD_NCH
#define DichLogisticBMD_NCH
#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>


#include "log_likelihoods.h"
#include "binomModels.h"



//void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd, void*)> math_func)
struct logistic_inequality {
	double BMD;
	double BMR;
	bool geq;
	double inequality;
};

double LOGISTIC_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);
double LOGISTIC_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);

class dich_logisticModelNC : public binomialBMD {

public:
	// hill dose response function
	dich_logisticModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int T) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 2);
		Eigen::MatrixXd one(temp.rows(), 1);
		newX << one,  temp;
		X = newX;
	};
	int    nParms() { return 2; };

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = LOGISTIC_A(theta(0, 0));
		double b = LOGISTIC_B(theta(1, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = LOGISTIC_MEAN(g, b, X(i, 1));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR, double isExtra) {
		double a = LOGISTIC_A(theta(0, 0));

		double rV = (isExtra) ? -1.0 : -1.0;

		//if gradient is specified we calculate the gradient
		if (grad) {
			if (isExtra) {
				grad[0] = 0.0;
			}
			else {
				grad[0] = 0.0;
			}
		}
		return rV;
	}

	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR) {
		double a = LOGISTIC_A(theta(0, 0));
		double b = LOGISTIC_B(theta(1, 0));
		double BMD = LOGISTIC_ADDED_RISK(a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double a   = LOGISTIC_A(theta(0, 0));
		double b   = LOGISTIC_B(theta(1, 0));
		double BMD = LOGISTIC_EXTRA_RISK(a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }


	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		double a = LOGISTIC_A(theta(0, 0));

		double Z, temp;
		if (isExtra) {
			Z = LOGISTIC_EXTRA_Z(a, BMR);
			//temp = (lb*log(BMD) - Z) *1.2; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);

		}else {
			Z = LOGISTIC_ADDED_Z(a, BMR);
			//temp = (lb*log(BMD/) - Z) * 1.2; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);

		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		double a = LOGISTIC_A(theta(0, 0));

		double Z, temp;
		if (isExtra) {
			Z = LOGISTIC_EXTRA_Z(a,  BMR);
			//temp = (Z - ub * log(BMD))*1.01; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			//x[1] = theta(1, 0) - temp;
		}
		else {
			Z = LOGISTIC_ADDED_Z(a,  BMR);
			//temp = (Z - ub * log(BMD))*1.01; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			//x[1] = theta(1, 0) - temp;
		}
		return x;

	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double a = LOGISTIC_A(theta(0, 0));
		double Z = LOGISTIC_ADDED_Z(a, BMR);

		double BETA = Z/BMD;

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double a = LOGISTIC_A(theta(0, 0));
		double Z = LOGISTIC_EXTRA_Z(a,  BMR);
		double BETA = Z/BMD;

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd XgivenD(double d) {
		Eigen::MatrixXd rV(1, 2); rV << 1.0, d; // by default the
													 // dose is the last parameter in X
		return rV;
	}

	double compute_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double *grad) {

		logistic_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, LOGISTIC_BMD_EXTRA_NC_INEQUALITY);
		}
		return LOGISTIC_BMD_EXTRA_NC_INEQUALITY(theta, &M);

	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		logistic_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, LOGISTIC_BMD_ADDED_NC_INEQUALITY);
		}
		return LOGISTIC_BMD_ADDED_NC_INEQUALITY(theta, &M);
	}

	int parameterRemoved() {
		return 1; // the Beta is the removed parameter for the LOGISTIC
	}
};

//Make sure the # defines are only for this file!



#endif
