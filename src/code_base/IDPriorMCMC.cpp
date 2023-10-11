#include "stdafx.h" // Precompiled header - does nothing if building R version

#ifdef R_COMPILATION
// necessary things to run in R
#include <RcppEigen.h>
#include <RcppGSL.h>
#else
#include <Eigen/Dense>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include "IDPriorMCMC.h"
///////////////////////////////////////////////////////////
// Function: neg_log_prior(Eigen::MatrixXd theta)
// Return the negative log prior for the prior specification
// Input: theta a kx1 matrix
// Output: Returns a real value based upon the prior for the given
//         paramer defined in prior_spec. Here prior_spec is a kx3
//         matrix.  The first column defines the type of prior 1-normal
//         2-log-normal, 0 improper.  The second colunn defines the mean.
//         For the improper prior, this is not used but it still defines an \
	//         option for the initialization.  The third column specifies the
//		   dispersion parameter.
double IDPriorMCMC::neg_log_prior(Eigen::MatrixXd theta, bool bound_check = false)
{
	double returnV = double(theta.rows()) * log(0.5 * M_2_SQRTPI * M_SQRT1_2);
	double mean = 0;
	double sd = 0;
	// loop over the prior specification in prior_spec
	// when  it is 1 - Normal Prior
	// when  it is 2 - Log normal prior.
	for (int i = 0; i < theta.rows(); i++)
	{
		int t = int(prior_spec(i, 0));
		if (theta(i, 0) < prior_spec(i, 3) ||
			theta(i, 0) > prior_spec(i, 4))
		{
			returnV = std::numeric_limits<double>::infinity();
			break;
		}
		switch (t)
		{
		case 1:
			mean = (theta(i, 0) - prior_spec(i, 1));
			sd = prior_spec(i, 2);
			returnV += -log(sd) - 0.5 * mean * mean / (sd * sd);

			break;
		case 2:
			mean = (log(theta(i, 0)) - prior_spec(i, 1));
			sd = prior_spec(i, 2);
			returnV += -log(sd) - log(theta(i, 0)) - 0.5 * mean * mean / (sd * sd);
			break;
		case 3:
			returnV -= log(0.5 * M_2_SQRTPI * M_SQRT1_2);
			returnV += log(gsl_ran_cauchy_pdf(theta(i,0) - prior_spec(i, 1), prior_spec(i,2)));
			break;
		case 4:
			returnV -= log(0.5 * M_2_SQRTPI * M_SQRT1_2);
			returnV += log(gsl_ran_gamma_pdf(theta(i, 0), pow(prior_spec(i,1)/prior_spec(i,2),2), pow(prior_spec(i ,2),2) / prior_spec(i, 1)));
			break;
		case 5:
		{
			returnV -= log(0.5 * M_2_SQRTPI * M_SQRT1_2);
			double gam = ((prior_spec(i,1) - prior_spec(i,3)) * (prior_spec(i,4) - prior_spec(i,1)))/(pow(prior_spec(i,2), 2)) - 3.0;
			double b = ((gam + 2.0) * prior_spec(i,1) - prior_spec(i,3) - prior_spec(i,4)) / gam;
			double alpha = 1.0 + gam * (b - prior_spec(i,3)) / (prior_spec(i,4) - prior_spec(i,3));
			double beta = 1.0 + gam * (prior_spec(i,4) - b) / (prior_spec(i,4) - prior_spec(i,3));
			returnV += (alpha - 1.0) * log(theta(i, 0) - prior_spec(i,3));
			returnV += (beta - 1.0) * log(prior_spec(i,4) - theta(i,0));
			returnV -= gsl_sf_lnbeta(alpha, beta);
			returnV -= (alpha + beta - 1.0) * log(prior_spec(i,4) - prior_spec(i,3));
		}
			break;
		default: // in the default case we remove all prior info
			returnV -= log(0.5 * M_2_SQRTPI * M_SQRT1_2);
			break;
		}
	}
	return -1.0 * returnV;
};

///////////////////////////////////////////////////////////
// Function: neg_log_prior(Eigen::MatrixXd theta)
// Return the negative log prior for the prior specification
// Input: theta a kx1 matrix
// Output: Returns a real value based upon the prior for the given
//         paramer defined in prior_spec. Here prior_spec is a kx3
//         matrix.  The first column defines the type of prior 1-normal
//         2-log-normal, 0 improper.  The second colunn defines the mean.
//         For the improper prior, this is not used but it still defines an \
//         option for the initialization.  The third column specifies the
//		   dispersion parameter.
Eigen::MatrixXd IDPriorMCMC::log_prior(Eigen::MatrixXd theta)
{
	double pi_const = log(0.5 * M_2_SQRTPI * M_SQRT1_2);
	Eigen::MatrixXd returnV(theta.rows(), 1);
	returnV.setZero();
	double mean = 0;
	double sd = 0;
	double mu = 0;
	// loop over the prior specification in prior_spec
	// when  it is 1 - Normal Prior
	// when  it is 2 - Log normal prior.
	for (int i = 0; i < theta.rows(); i++)
	{
		int t = int(prior_spec(i, 0));
		switch (t)
		{
		case 1:
			mean = (theta(i, 0) - prior_spec(i, 1));
			sd = prior_spec(i, 2);
			// returnV(i,0) = // - 0.5*mean*mean / (sd*sd);// +pi_const -log(sd);
			returnV(i, 0) = -1 / (sd * sd);
			break;
		case 2:
			mean = (log(theta(i, 0)) - prior_spec(i, 1));
			sd = prior_spec(i, 2);
			mu = prior_spec(i, 1);
			returnV(i, 0) = (exp(sd * sd) - 1);
			returnV(i, 0) *= exp(2 * mu + sd * sd);
			returnV(i, 0) = -1 / returnV(i, 0);
			// returnV(i,0) = - 0.5*mean*mean / (sd*sd) ;// -log(sd) - log(theta(i, 0)) + pi_const;
			break;
		case 3: case 4: case 5:
			sd = prior_spec(i,2);
			returnV(i,0) = -1.0 / (sd * sd);
			break;
		default: // in the default case we remove all prior info
			returnV(i, 0) = fabs(0.0);
			break;
		}
	}
	return -1 * returnV.asDiagonal();
};

//////////////////////////////////////////////////////
// Function: prior_mean()
// Purpose returns the prior_means
// Output: Returns the prior mean for normal and median for the log normal distributions
//		   for improper prior it returns the value that is in column 2, which
//         can be used for initial starting points.
Eigen::MatrixXd IDPriorMCMC::prior_mean()
{
	Eigen::MatrixXd pmean(prior_spec.rows(), 1);
	pmean.setZero(); 

	for (int i = 0; i < prior_spec.rows(); i++)
	{
		int t = int(prior_spec(i, 0));
		switch (t)
		{
		case 2:
			pmean(i, 0) = exp(prior_spec(i, 1));
			break;
		default: // in the default case we remove all prior info
			pmean(i, 0) = prior_spec(i, 1);
			break;
		}
	}
	return pmean;
};
