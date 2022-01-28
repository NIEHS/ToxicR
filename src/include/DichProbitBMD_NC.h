//File: DichPROBITlBMD_NC.h
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
#define PROBIT_A(X) X
#define PROBIT_B(X) X
#define ADDED_Q(A,BMR) (BMR+gsl_cdf_gaussian_P(A,1.0))
#define EXTRA_Q(A,BMR) (BMR*(1.0-gsl_cdf_gaussian_P(A,1.0))+gsl_cdf_gaussian_P(A,1.0))
#define PROBIT_ADDED_Z(A,BMR) (gsl_cdf_gaussian_Pinv(ADDED_Q(A,BMR),1.0) - A)
#define PROBIT_EXTRA_Z(A,BMR) (gsl_cdf_gaussian_Pinv(EXTRA_Q(A,BMR),1.0) - A)
#define PROBIT_ADDED_RISK(A,B,BMR) PROBIT_ADDED_Z(A,BMR)/B
#define PROBIT_EXTRA_RISK(A,B,BMR) PROBIT_EXTRA_Z(A,BMR)/B
#define PROBIT_MEAN(A,B,D) gsl_cdf_gaussian_P(A+B*D,1.0)



#ifndef DichProbitBMD_NCHB
#define DichProbitBMD_NCHB


#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include "log_likelihoods.h"
#include "binomModels.h"


//void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd, void*)> math_func)
struct probit_inequality {
	double BMD;
	double BMR;
	bool geq;
	double inequality;
};

double PROBIT_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);
double PROBIT_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);

class dich_probitModelNC : public binomialBMD {

public:
	// hill dose response function
	dich_probitModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int Deg) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 2);
		Eigen::MatrixXd one(temp.rows(), 1) ;
		newX << one, temp;
		X = newX;
	};
	int    nParms() { return 2; };

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = PROBIT_A(theta(0, 0));
		double b = PROBIT_B(theta(1, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = PROBIT_MEAN(g, b, X(i, 1));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR, double isExtra) {
		double a = PROBIT_A(theta(0, 0));

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
		double a = PROBIT_A(theta(0, 0));
		double b = PROBIT_B(theta(1, 0));
		double BMD = PROBIT_ADDED_RISK(a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double a   = PROBIT_A(theta(0, 0));
		double b   = PROBIT_B(theta(1, 0));
		double BMD = PROBIT_EXTRA_RISK(a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }


	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		double a = PROBIT_A(theta(0, 0));

		double Z, temp;
		if (isExtra) {
			Z = PROBIT_EXTRA_Z(a, BMR);
			//temp = (lb*log(BMD) - Z) *1.2; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);

		}else {
			Z = PROBIT_ADDED_Z(a, BMR);
			//temp = (lb*log(BMD/) - Z) * 1.2; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);

		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		double a = PROBIT_A(theta(0, 0));

		double Z, temp;
		if (isExtra) {
			Z = PROBIT_EXTRA_Z(a,  BMR);
			//temp = (Z - ub * log(BMD))*1.01; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			//x[1] = theta(1, 0) - temp;
		}
		else {
			Z = PROBIT_ADDED_Z(a,  BMR);
			//temp = (Z - ub * log(BMD))*1.01; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			//x[1] = theta(1, 0) - temp;
		}
		return x;

	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double a = PROBIT_A(theta(0, 0));
		double Z = PROBIT_ADDED_Z(a, BMR);

		double BETA = Z/BMD;

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double a = PROBIT_A(theta(0, 0));
		double Z = PROBIT_EXTRA_Z(a,  BMR);
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

		probit_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, PROBIT_BMD_EXTRA_NC_INEQUALITY);
		}
		return PROBIT_BMD_EXTRA_NC_INEQUALITY(theta, &M);

	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		probit_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, PROBIT_BMD_ADDED_NC_INEQUALITY);
		}
		return PROBIT_BMD_ADDED_NC_INEQUALITY(theta, &M);
	}

	int parameterRemoved() {
		return 1; // the Beta is the removed parameter for the PROBIT
	}
};




#endif
