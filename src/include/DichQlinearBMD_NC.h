//File: DichQLINEARlBMD_NC.h
//Purpose: Creates a Dichotomous QuantalLinear BMD Model
//         across different statisticalmodels. Log_likelihood classes
//         define specific probability models.
//Creator: Matt Wheeler
//Date   : 1/25/2018
//Changes:
//
//
//
#define QLINEAR_G(X) 1.0/(1.0 + exp(-X))
#define QLINEAR_B(X) X
#define QLINEAR_EXTRA_Z(G,BMR) -log(1.0-BMR)
#define QLINEAR_ADDED_Z(G,BMR) -log(1.0-BMR/(1.0-G))
#define QLINEAR_ADDED_RISK(G,B,BMR) QLINEAR_ADDED_Z(G,BMR)/B
#define QLINEAR_EXTRA_RISK(G,B,BMR) QLINEAR_EXTRA_Z(G,BMR)/B
#define QLINEAR_MEAN(G,B,D) ((D) <= (0.0)) ? (G) : (G+(1.0-G)*(1.0-exp(-B*D)))


#ifndef DichQlinearBMD_NCH
#define DichQlinearBMD_NCH
#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif

#include "log_likelihoods.h"
#include "binomModels.h"

double QLINEAR_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);
double QLINEAR_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data);

//void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd, void*)> math_func)
struct log_qlinear_inequality {
	double BMD;
	double BMR;
	bool geq;
	double inequality;
};



class dich_qlinearModelNC : public binomialBMD {

public:

	dich_qlinearModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int Deg) : binomialBMD(tY, tX) {
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 2);
		Eigen::MatrixXd one(temp.rows(), 1);
		newX << one, temp;
		X = newX;
	};
	
	int    nParms() { return 2; };
	
	virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D){
	     Eigen::MatrixXd returnV(D.rows(), 2);
	     Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1);
	     returnV << one , D.col(0); 
	      
	     return returnV; 
	     
	}
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {

		double g = QLINEAR_G(theta(0, 0));
		double b = QLINEAR_B(theta(1, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = QLINEAR_MEAN(g, b, X(i, 1));

		return   p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double * grad, double BMR, double isExtra) {
		double g = QLINEAR_G(theta(0, 0));

		double rV = (isExtra) ? -1.0 : BMR / (1 - g) - 1;

		//if gradient is specified we calculate the gradient
		if (grad) {
			if (isExtra) {
				grad[0] = 0.0;
			}
			else {
				grad[0] = -BMR * exp(theta(0, 0)) / pow(BMR + exp(theta(0, 0)), 2.0);
			}
		}
		return rV;
	}

	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR) {
		double g = QLINEAR_G(theta(0, 0));
		double b = QLINEAR_B(theta(1, 0));
		double BMD = QLINEAR_ADDED_RISK(g, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR) {
		double g = QLINEAR_G(theta(0, 0));
		double b = QLINEAR_B(theta(1, 0));
		double BMD = QLINEAR_EXTRA_RISK(g, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }


	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		std::vector<double> x(theta.rows());
		double g = QLINEAR_G(theta(0, 0));

		double Z, temp;
		if (isExtra) {
			Z = QLINEAR_EXTRA_Z(g, BMR);
			//temp = (lb*log(BMD) - Z) *1.2; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);

		}else {
			Z = QLINEAR_ADDED_Z(g, BMR);
			//temp = (lb*log(BMD/) - Z) * 1.2; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);

		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {

		std::vector<double> x(theta.rows());
		double g = QLINEAR_G(theta(0, 0));

		double Z, temp;
		if (isExtra) {
			Z = QLINEAR_EXTRA_Z(g,  BMR);
			//temp = (Z - ub * log(BMD))*1.01; // add a little extra add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			//x[1] = theta(1, 0) - temp;
		}
		else {
			Z = QLINEAR_ADDED_Z(g,  BMR);
			//temp = (Z - ub * log(BMD))*1.01; // add a little extra so the constraint is satisfied
			x[0] = theta(0, 0);
			//x[1] = theta(1, 0) - temp;
		}
		return x;

	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double g = QLINEAR_G(theta(0, 0));
		double Z = QLINEAR_ADDED_Z(g, BMR);

		double BETA = Z/BMD;

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0); newM(1, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR) {
		double g = QLINEAR_G(theta(0, 0));
		double Z = QLINEAR_EXTRA_Z(g,  BMR);
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

		log_qlinear_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, QLINEAR_BMD_EXTRA_NC_INEQUALITY);
		}
		return QLINEAR_BMD_EXTRA_NC_INEQUALITY(theta, &M);

	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double * grad) {
		log_qlinear_inequality M;
		M.inequality = inequality; M.BMD = BMD;
		M.BMR = BMR; M.geq = geq;

		if (grad) {
			gradient(theta, grad, &M, QLINEAR_BMD_ADDED_NC_INEQUALITY);
		}
		return QLINEAR_BMD_ADDED_NC_INEQUALITY(theta, &M);
	}

	int parameterRemoved() {
		return 1; // the Beta is the removed parameter for the QLINEAR
	}
};

//Make sure the # defines are only for this file!
#endif
