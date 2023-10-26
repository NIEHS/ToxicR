// File: DichGAMMAlBMD_NC.h
// Purpose: Creates a Dichotomous Hill BMD Model
//          across different statisticalmodels. Log_likelihood classes
//          define specific probability models.
// Creator: Matt Wheeler
// Date   : 1/5/2018
// Changes:
//
//
//

#define GAMMA_G(X) 1.0 / (1.0 + exp(-X))
#define GAMMA_A(X) X
#define GAMMA_B(X) X
#define GAMMA_ADDED_Z(G, A, BMR) gsl_cdf_gamma_Pinv(BMR / (1 - G), A, 1)
#define GAMMA_EXTRA_Z(G, A, BMR) gsl_cdf_gamma_Pinv(BMR, A, 1)
#define GAMMA_ADDED_RISK(G, A, B, BMR) GAMMA_ADDED_Z(G, A, BMR) / B
#define GAMMA_EXTRA_RISK(G, A, B, BMR) GAMMA_EXTRA_Z(G, A, BMR) / B
#define GAMMA_MEAN(G, A, B, D) ((D) <= (0.0)) ? (G) : (G + (1 - G) * (gsl_cdf_gamma_P(B * D, A, 1)))

#pragma once
#ifndef DichGammaBMD_NCH
#define DichGammaBMD_NCH

#ifdef R_COMPILATION
// necessary things to run in R
#include <RcppGSL.h>
#include <RcppEigen.h>
#else
#include <Eigen/Dense>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#include "log_likelihoods.h"
#include "binomModels.h"


class ErrorHandler {
	public: 

	static void gsl_err_gamma_extra_z_handler(const char * reason,
	const char * file, 
	int line, 
	int gsl_errno
	) {
		std::cerr << "GSL error in file " << file << ", line " << line << ", reason " << reason << std::endl;
		if (gsl_errno == GSL_EDOM)  {
			std::cerr << "Domain Error!" << std::endl;
		} else if (gsl_errno == GSL_EINVAL) {
			std::cerr << "Invalid argument!" << std::endl;
		} else if (gsl_errno == GSL_EROUND) {
			std::cerr << "Round off error!" << std::endl;
		}
	}
};

// void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd, void*)> math_func)
struct log_gamma_inequality
{
	double BMD;
	double BMR;
	bool geq;
	double inequality;
};

double GAMMA_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void *data);
double GAMMA_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void *data);

class dich_gammaModelNC : public binomialBMD
{

public:
	// hill dose response function
	dich_gammaModelNC(Eigen::MatrixXd tY, Eigen::MatrixXd tX, int Deg) : binomialBMD(tY, tX)
	{
		Eigen::MatrixXd temp = X;
		Eigen::MatrixXd newX(temp.rows(), 3);
		Eigen::MatrixXd one(temp.rows(), 1);
		one.setZero(); 
		newX.setZero(); 
		newX << one, one, temp;
		X = newX;
	};
	int nParms() { return X.cols(); };
	virtual Eigen::MatrixXd convertDataMatrix(Eigen::MatrixXd D)
	{
		Eigen::MatrixXd returnV(D.rows(), 3);
		Eigen::MatrixXd one = Eigen::MatrixXd::Ones(D.rows(), 1);
		returnV << one, one, D.col(0);
		return returnV;
	}
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta)
	{
		return mean(theta, X);
	};

	Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X)
	{

		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		double b = GAMMA_B(theta(2, 0));

		Eigen::MatrixXd p(X.rows(), 1);

		for (int i = 0; i < X.rows(); i++)
			p(i, 0) = GAMMA_MEAN(g, a, b, X(i, 2));

		return p;
	};

	virtual double BMR_CONSTRAINT(Eigen::MatrixXd theta, double *grad, double BMR, double isExtra)
	{
		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		double b = GAMMA_B(theta(2, 0));

		double rV = (isExtra) ? -1.0 : BMR / (1 - g) - 1;

		// if gradient is specified we calculate the gradient
		if (grad)
		{
			if (isExtra)
			{
				grad[0] = 0.0;
				grad[1] = 0.0;
			}
			else
			{
				grad[0] = -BMR * exp(theta(0, 0)) / pow(BMR + exp(theta(0, 0)), 2.0);
				grad[1] = 0.0;
			}
		}
		return rV;
	}

	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_ADDED_NC(Eigen::MatrixXd theta, double BMR)
	{
		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		double b = GAMMA_B(theta(2, 0));

		double BMD = GAMMA_ADDED_RISK(g, a, b, BMR);
		return BMD;
	}
	// Computes the BMD with a given BMR assuming no covariates NC option
	double compute_BMD_EXTRA_NC(Eigen::MatrixXd theta, double BMR)
	{
		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		double b = GAMMA_B(theta(2, 0));

		double BMD = GAMMA_EXTRA_RISK(g, a, b, BMR);
		return BMD;
	}

	// returns true if you can turn the BMD as a function of the
	// other parameters false if it is a equality constraint
	// defaults to false
	bool fixedNC_BMDformula() { return true; }

	std::vector<double> fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb)
	{
		std::vector<double> x(theta.rows());
		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		double temp_a = a;
		double upper, lower;
		double Z, temp;
		if (isExtra)
		{
			Z = GAMMA_EXTRA_Z(g, a, BMR);
			temp = Z / BMD + fabs(lb - Z / BMD) * 1.01; // define a beta that is "within" the constraints
			x[0] = theta(0, 0);

			while (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD < 0)
			{
				temp_a *= .95;
			}

			upper = temp_a;
			temp_a = a;
			while (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD > 0)
			{
				temp_a *= 1.05;
			}
			lower = temp_a;
			temp_a = 0.5 * (upper + lower);

			while (fabs(temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD) > 1e-2)
			{
				if (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD < 0)
				{
					lower = temp_a;
				}
				else
				{
					upper = temp_a;
				}
				temp_a = 0.5 * (upper + lower);
			}
			x[1] = temp_a;
		}
		else
		{
			Z = GAMMA_ADDED_Z(g, a, BMR);
			temp = Z / BMD + fabs(lb - Z / BMD) * 1.01; // define a beta that is "within" the constraints
			x[0] = theta(0, 0);

			while (temp - GAMMA_ADDED_Z(g, temp_a, BMR) / BMD < 0)
			{
				temp_a *= .95;
			}

			upper = temp_a;
			temp_a = a;
			while (temp - GAMMA_ADDED_Z(g, temp_a, BMR) / BMD > 0)
			{
				temp_a *= 1.05;
			}
			lower = temp_a;
			temp_a = 0.5 * (upper + lower);

			while (fabs(temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD) > 1e-2)
			{
				if (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD < 0)
				{
					lower = temp_a;
				}
				else
				{
					upper = temp_a;
				}
				temp_a = 0.5 * (upper + lower);
			}
			x[1] = temp_a;
		}
		return x;
	}

	std::vector<double> fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub)
	{

		std::vector<double> x(theta.rows());
		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		double temp_a = a;

		double upper, lower;
		double Z, temp;
		if (isExtra)
		{
			Z = GAMMA_EXTRA_Z(g, a, BMR);
			temp = Z / BMD - fabs(ub - Z / BMD) * 1.01; // define a beta that is "within" the constraints
			x[0] = theta(0, 0);
			// set up a zero in procedure to get a good alpha estimate

			while (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD < 0)
			{
				temp_a *= .95;
			}

			upper = temp_a;
			temp_a = a;
			while (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD > 0)
			{
				temp_a *= 1.05;
			}
			lower = temp_a;
			temp_a = 0.5 * (upper + lower);

			while (fabs(temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD) > 1e-2)
			{
				if (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD < 0)
				{
					lower = temp_a;
				}
				else
				{
					upper = temp_a;
				}
				temp_a = 0.5 * (upper + lower);
			}
			x[1] = temp_a;
		}
		else
		{
			Z = GAMMA_ADDED_Z(g, a, BMR);
			temp = Z / BMD - fabs(ub - Z / BMD) * 1.01; // define a beta that is "within" the constraints
			x[0] = theta(0, 0);
			// set up a zero in procedure to get a good alpha estimate

			while (temp - GAMMA_ADDED_Z(g, temp_a, BMR) / BMD < 0)
			{
				temp_a *= .95;
			}

			upper = temp_a;
			temp_a = a;
			while (temp - GAMMA_ADDED_Z(g, temp_a, BMR) / BMD > 0)
			{
				temp_a *= 1.05;
			}
			lower = temp_a;
			temp_a = 0.5 * (upper + lower);

			while (fabs(temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD) > 1e-2)
			{
				if (temp - GAMMA_EXTRA_Z(g, temp_a, BMR) / BMD < 0)
				{
					lower = temp_a;
				}
				else
				{
					upper = temp_a;
				}
				temp_a = 0.5 * (upper + lower);
			}
			x[1] = temp_a;
		}
		return x;
	}

	Eigen::MatrixXd beta_BMD_AddedNC(Eigen::MatrixXd theta, double BMD, double BMR)
	{
		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));

		double Z = GAMMA_ADDED_Z(g, a, BMR);
		double BETA = Z / BMD;

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0);
		newM(1, 0) = theta(1, 0);
		newM(2, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd beta_BMD_ExtraNC(Eigen::MatrixXd theta, double BMD, double BMR)
	{

		double g = GAMMA_G(theta(0, 0));
		double a = GAMMA_A(theta(1, 0));
		// clamp BMR to 0-1 (avoiding memory leaks)
		BMR = std::max(1e-6, std::min(BMR, 1.0 - 1e-6));
		
		// gsl_set_error_handler(ErrorHandler::gsl_err_gamma_extra_z_handler);

		double Z = GAMMA_EXTRA_Z(g, a, BMR);
		

		double BETA = Z / BMD;

		Eigen::MatrixXd newM(theta.rows() + 1, 1);
		newM(0, 0) = theta(0, 0);
		newM(1, 0) = theta(1, 0);
		newM(2, 0) = BETA;
		return newM;
	}

	Eigen::MatrixXd XgivenD(double d)
	{
		Eigen::MatrixXd rV(1, 3);
		rV << 1.0, 1.0, d; // by default the
						   // dose is the last parameter in X
		return rV;
	}

	double compute_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double *grad)
	{

		log_gamma_inequality M;
		M.inequality = inequality;
		M.BMD = BMD;
		M.BMR = BMR;
		M.geq = geq;

		if (grad)
		{
			gradient(theta, grad, &M, GAMMA_BMD_EXTRA_NC_INEQUALITY);
		}
		return GAMMA_BMD_EXTRA_NC_INEQUALITY(theta, &M);
	}

	double compute_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq, double *grad)
	{
		log_gamma_inequality M;
		M.inequality = inequality;
		M.BMD = BMD;
		M.BMR = BMR;
		M.geq = geq;

		if (grad)
		{
			gradient(theta, grad, &M, GAMMA_BMD_ADDED_NC_INEQUALITY);
		}
		return GAMMA_BMD_ADDED_NC_INEQUALITY(theta, &M);
	}

	int parameterRemoved()
	{
		return 2; // the Beta is the removed parameter for the GAMMA
	}
};

// Make sure the # defines are only for this file!
#endif
