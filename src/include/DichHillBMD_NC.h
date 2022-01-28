//File: DichHillBMD_NC.h
//Purpose: Creates a Dichotomous Hill BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/25/2018
//Changes:
//  11/09/2018, Matt Wheeler
//   Fixed the derivatives for the added risk inequality, and thus
//   the MMA optimization.  This solves the "slow runtime" issue.
//
//
//
#pragma once
#ifndef DichHillBMD_NCH
#define DichHillBMD_NCH

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif


#include "log_likelihoods.h"
#include "binomModels.h"

// The background and asymptote terms are transformed to a logisitc distribution to
// improve optimization and need to be re-transformed before being used.
#define HILL_G(X) 1.0/(1.0 + exp(-X))
#define HILL_N(X) 1.0/(1.0 + exp(-X))
#define HILL_A(X) X
#define HILL_B(X) X
#define HILL_ADDED_Z(G,N,A,B,BMR) -A - log((1 - G)*N / BMR - 1.0)
#define HILL_EXTRA_Z(G,N,A,B,BMR)  -A - log(N / BMR - 1.0);
#define HILL_ADDED_RISK(G,N,A,B,BMR) exp((-A - log((1 - G)*N / BMR - 1.0))/B)
#define HILL_EXTRA_RISK(G,N,A,B,BMR) exp((-A - log(N / BMR - 1.0))/B)
// To improve optimization, the model function was intentionally re-parameterized
// such that the background = G, not G*V
#define HILL_MEAN(G,N,A,B,D) ((D) <= (0.0)) ? (G) : (G+(1-G)*N/(1+exp(-A-B*log(D))))


class dich_hillModelNC : public binomialBMD {

public:
	// hill dose response function
	dich_hillModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int T) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 3);
		Eigen::MatrixXd one(temp.rows(), 1);
		newX << one, one, temp;
		X = newX;
	};
	int    nParms() { return X.cols() + 1; }; // assumes n x 3 matrix +1 is four
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};
	
	virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
	     Eigen::MatrixXd returnV(D.rows(), 3);
	     Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1) ;
	     returnV << one , one , D.col(0); 
	     return returnV; 
	     
	}
	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {
		std::vector<double> x(theta.rows());
		double g = HILL_G(theta(0, 0)); // + 0.5;
		double n = HILL_N(theta(1, 0));// + 0.5;
		double a = HILL_A(theta(2, 0));
		double Z, temp;
		if (isExtra) {
			Z = HILL_EXTRA_Z(g, n, a, NULL, BMR);
			temp = (Z - ub * log(BMD))*1.1; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0);
			x[2] = theta(2, 0) - temp;
		}
		else {
			Z = HILL_ADDED_Z(g,n, a, NULL, BMR);
			temp = temp = (Z - ub * log(BMD))*1.1; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0);
			x[2] = theta(2, 0) - temp;
		}
		return x;
	}

	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		double g = HILL_G(theta(0, 0)); // + 0.5;
		double n = HILL_N(theta(1, 0));// + 0.5;
		double a = HILL_A(theta(2, 0));
		double Z, temp;
		if (isExtra) {
			Z = HILL_EXTRA_Z(g,n,a,NULL, BMR);
			temp = (lb*log(BMD) - Z) * 2; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0);
			x[2] = theta(2, 0) - temp;
		}
		else {
			Z = HILL_ADDED_Z(g, n, a,NULL, BMR);
			temp = (lb*log(BMD) - Z) * 2; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			x[1] = theta(1, 0);
			x[2] = theta(2, 0) - temp;
		}
		return x;
	}

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
   
		double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2,0)); double b = HILL_B(theta(3, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = HILL_MEAN(g, n, a, b, X(i, 2));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR , double isExtra) {
		double g = HILL_G(theta(0, 0)); // + 0.5;
		double n = HILL_N(theta(1, 0));// + 0.5;
		//double a = HILL_A(theta(2, 0)); //double b = HILL_B(theta(3, 0));

		double rV = (isExtra)? -n/BMR + 1.0 : -((1.0 - g)*n) / BMR + 1.0;

		//if gradient is specified we calculate the gradient
		if (grad) {
			if (isExtra) {
				grad[0] = 0.0; grad[1] = -1.0 / BMR; grad[2] = 0.0;
			}
			else {
				grad[0] = n / BMR; grad[1] = -(1.0 - g) / BMR; grad[2] = 0.0;
			}
		}
		return rV;
	}

	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR) {
		double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2, 0)); double b = HILL_B(theta(3, 0));

		double BMD = HILL_ADDED_RISK(g, n, a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
	//	double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2, 0)); double b = HILL_B(theta(3, 0));

		double BMD = HILL_EXTRA_RISK(g, n, a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2, 0));

		double Z = HILL_ADDED_Z(g, n, a, b, BMR);
		double BETA = Z / log(BMD);

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = theta(1, 0);
		newM(2, 0) = theta(2, 0);  newM(3, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {

	//	double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2, 0));

		double Z = HILL_EXTRA_Z(g, n, a, b, BMR);

		double BETA = Z / log(BMD);

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = theta(1, 0);
		newM(2, 0) = theta(2, 0);  newM(3, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd XgivenD(double d) {
		Eigen::MatrixXd rV(1, 3); rV << 1.0, 1.0, d; // by default the
													 // dose is the last parameter in X
		return rV;
	}

	double compute_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double *grad) {

	//	double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2, 0));
		double Z = HILL_EXTRA_Z(g, n, a, b, BMR);

		//grad is non null compute the gradient of the inequality
		// constraint
		if (grad) {
			grad[0] = 0.0; //gamma
			grad[1] = -1.0*(n / BMR)*exp(theta(1, 0)) / ((exp(theta(1, 0)) + 1)*((n / BMR - 1.0)*exp(theta(1, 0)) - 1));
			grad[2] = -1.0; // alpha
		}

		Z = Z / log(BMD);

		double rV = 0.0;
		if (geq) { // greater than or equal
			rV = inequality - Z + 1e-6;
			if (grad) {
				grad[1] = grad[1] * (-1.0 / log(BMD));
				grad[2] = grad[2] * (-1.0 / log(BMD));
			}
		}
		else {
			rV = Z - inequality + 1e-6;
			if (grad) {
				grad[1] = grad[1] * (-1.0 / log(BMD));
				grad[2] = grad[2] * (-1.0 / log(BMD));
			}

		}

		return rV;
	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		double g = HILL_G(theta(0, 0));
		double n = HILL_N(theta(1, 0));
		double a = HILL_A(theta(2, 0));

		double Z = HILL_ADDED_Z(g, n, a, b, BMR);
		Z = Z / log(BMD);
		double rV = 0.0;

		//grad is non null compute the gradient
		if (grad) {
			grad[0] = -(n/(n*g-n+BMR))*(exp(theta(1,0))/((exp(theta(1,0)+1)*(exp(theta(1,0)+1))))); //gamma
			grad[1] =  (g-1)/((g-1)*n+BMR)*(exp(theta(1,0))/((exp(theta(1,0)+1)*(exp(theta(1,0)+1)))));
			grad[2] = -1.0; // alpha
		}

		if (geq) { // greater than or equal
			rV = inequality - Z ;
			if (grad) {
				grad[0] *= (-1.0 / log(BMD));
				grad[1] *= (-1.0 / log(BMD));
				grad[2] *= (-1.0 / log(BMD));
			}
		}
		else {
			rV = Z - inequality ;
			if (grad) {
				grad[0] *= (1.0 / log(BMD));
				grad[1] *= (1.0 / log(BMD));
				grad[2] *= (1.0 / log(BMD));
			}

		}

		return rV;
	}

	int parameterRemoved() {
		return 3; // the Beta is the removed parameter for the Hill
	}
};

//Make sure the # defines are only for this file!
#undef HILL_G
#undef HILL_N
#undef HILL_A
#undef HILL_B
#undef HILL_EXTRA_Z
#undef HILL_ADDED_RISK
#undef HILL_EXTRA_RISK
#endif
