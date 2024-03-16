/*
 * Copyright 2020  US HHS, NIEHS 
 * Author Matt Wheeler 
 * e-mail: <matt.wheeler@nih.gov> 
 *
 *Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 *and associated documentation files (the "Software"), to deal in the Software without restriction,
 *including without limitation the rights to use, copy, modify, merge, publish, distribute,
 *sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 *is furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all copies
 *or substantial portions of the Software.

 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 *CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 *OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 */
#ifndef bmd_calculateH
#define bmd_calculateH

//#include "stdafx.h"
#include <chrono>
#include <cmath>
#ifndef WIN32
#include <cfloat>
#endif
#include <chrono>
#include <vector>
#include <iostream>
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
#include <gsl/gsl_cdf.h>

#include "gradient.h"
#include "normal_EXP_NC.h"
#include "dBMDstatmod.h"
#include "cBMDstatmod.h"
#include "bmdStruct.h"
#include "continuous_clean_aux.h"

class bmd_cdf
{

public:
	bmd_cdf() : probs(), BMD()
	{
		spline_bmd_cdf = NULL;
		spline_bmd_inv = NULL;
		acc_bmd_cdf = NULL;
		acc_bmd_inv = NULL;
		multiple = 1.0;
		max_BMD = 0.0;
		min_BMD = 0.0;

		min_prob = 0.0;
		max_prob = 0.0;
	};

	bmd_cdf &operator=(const bmd_cdf &M)
	{

		probs = M.probs;
		BMD = M.BMD;
		multiple = M.multiple;
		max_BMD = M.max_BMD;
		min_BMD = M.min_BMD;
		max_prob = M.max_prob;
		min_prob = M.min_prob;
		int error = 0;

		if (probs.size() == BMD.size() && BMD.size() > 0)
		{
			acc_bmd_inv = gsl_interp_accel_alloc();
			acc_bmd_cdf = gsl_interp_accel_alloc();
			spline_bmd_inv = gsl_spline_alloc(gsl_interp_steffen, BMD.size());
			spline_bmd_cdf = gsl_spline_alloc(gsl_interp_steffen, BMD.size());
			error = gsl_spline_init(spline_bmd_inv, probs.data(), BMD.data(), BMD.size());

			if (error)
			{
				if (spline_bmd_inv)
				{
					gsl_spline_free(spline_bmd_inv);
				}
				if (spline_bmd_cdf)
				{
					gsl_spline_free(spline_bmd_cdf);
				}
				if (acc_bmd_cdf)
				{
					gsl_interp_accel_free(acc_bmd_cdf);
				}
				if (acc_bmd_inv)
				{
					gsl_interp_accel_free(acc_bmd_inv);
				}
				spline_bmd_inv = NULL;
				acc_bmd_inv = NULL;
				spline_bmd_inv = NULL;
				acc_bmd_inv = NULL;
			}
			else
			{
				error = gsl_spline_init(spline_bmd_cdf, BMD.data(), probs.data(), BMD.size());

				if (error)
				{
					if (spline_bmd_inv)
					{
						gsl_spline_free(spline_bmd_inv);
					}
					if (spline_bmd_cdf)
					{
						gsl_spline_free(spline_bmd_cdf);
					}
					if (acc_bmd_cdf)
					{
						gsl_interp_accel_free(acc_bmd_cdf);
					}
					if (acc_bmd_inv)
					{
						gsl_interp_accel_free(acc_bmd_inv);
					}
					// Rcpp::Rcout << "error: %s\n" << gsl_strerror (error) << endl;
					spline_bmd_cdf = NULL;
					acc_bmd_cdf = NULL;
					spline_bmd_inv = NULL;
					acc_bmd_inv = NULL;
				}
			}
		}
		return *this;
	}

	bmd_cdf(const bmd_cdf &M)
	{
		*this = M;
	}

	bmd_cdf(std::vector<double> tx, std::vector<double> ty) : probs(tx), BMD(ty)
	{
		int error;
		multiple = 1.0;
		max_prob = *max_element(probs.begin(), probs.end());
		min_prob = *min_element(probs.begin(), probs.end());

		max_BMD = *max_element(BMD.begin(), BMD.end());
		min_BMD = *min_element(BMD.begin(), BMD.end());

		if (probs.size() == BMD.size() && BMD.size() > 0)
		{
			acc_bmd_inv = gsl_interp_accel_alloc();
			acc_bmd_cdf = gsl_interp_accel_alloc();
			spline_bmd_inv = gsl_spline_alloc(gsl_interp_steffen, BMD.size());
			spline_bmd_cdf = gsl_spline_alloc(gsl_interp_steffen, BMD.size());
			error = gsl_spline_init(spline_bmd_inv, probs.data(), BMD.data(), BMD.size());

			if (error)
			{
				if (spline_bmd_inv)
				{
					gsl_spline_free(spline_bmd_inv);
				}
				if (spline_bmd_cdf)
				{
					gsl_spline_free(spline_bmd_cdf);
				}
				if (acc_bmd_cdf)
				{
					gsl_interp_accel_free(acc_bmd_cdf);
				}
				if (acc_bmd_inv)
				{
					gsl_interp_accel_free(acc_bmd_inv);
				}
				spline_bmd_inv = NULL;
				acc_bmd_inv = NULL;
				spline_bmd_inv = NULL;
				acc_bmd_inv = NULL;
			}
			else
			{
				error = gsl_spline_init(spline_bmd_cdf, BMD.data(), probs.data(), BMD.size());

				if (error)
				{
					if (spline_bmd_inv)
					{
						gsl_spline_free(spline_bmd_inv);
					}
					if (spline_bmd_cdf)
					{
						gsl_spline_free(spline_bmd_cdf);
					}
					if (acc_bmd_cdf)
					{
						gsl_interp_accel_free(acc_bmd_cdf);
					}
					if (acc_bmd_inv)
					{
						gsl_interp_accel_free(acc_bmd_inv);
					}
					// Rcpp::Rcout << "error: %s\n" << gsl_strerror (error) << endl;
					spline_bmd_cdf = NULL;
					acc_bmd_cdf = NULL;
					spline_bmd_inv = NULL;
					acc_bmd_inv = NULL;
				}
			}
		}
	}

	~bmd_cdf()
	{
		if (spline_bmd_inv)
		{
			gsl_spline_free(spline_bmd_inv);
		}
		if (spline_bmd_cdf)
		{
			gsl_spline_free(spline_bmd_cdf);
		}
		if (acc_bmd_cdf)
		{
			gsl_interp_accel_free(acc_bmd_cdf);
		}
		if (acc_bmd_inv)
		{
			gsl_interp_accel_free(acc_bmd_inv);
		}
		acc_bmd_inv = NULL;
		acc_bmd_cdf = NULL;
		spline_bmd_cdf = NULL;
		spline_bmd_inv = NULL;
	}

	double inv(double p)
	{

		if (spline_bmd_cdf != NULL && acc_bmd_cdf != NULL)
		{

			if (p > min_prob && p < max_prob)
			{
				return gsl_spline_eval(spline_bmd_inv, p, acc_bmd_inv) * multiple;
			}
			if (p < min_prob)
				return 0;
			if (p > max_prob)
				return std::numeric_limits<double>::infinity();
		}
		return NAN;
	}

	double P(double dose)
	{
		if (spline_bmd_cdf != NULL && acc_bmd_cdf != NULL)
		{
			if (dose / multiple > min_BMD && dose / multiple < max_BMD)
				return gsl_spline_eval(spline_bmd_cdf, dose / multiple, acc_bmd_cdf);
			if (dose / multiple < min_BMD)
				return 0.0;
			if (dose / multiple > max_BMD)
				return 1.0;
		}
		return NAN;
	}
	void set_multiple(double m)
	{
		if (m > 0.0 && std::isnormal(m))
		{
			multiple = m;
		}
		else
		{
			multiple = 1.0;
		}
	}

private:
	double min_BMD;
	double max_BMD;
	double multiple;
	double min_prob;
	double max_prob;
	std::vector<double> probs;
	std::vector<double> BMD;
	gsl_interp_accel *acc_bmd_cdf;
	gsl_spline *spline_bmd_cdf;
	gsl_interp_accel *acc_bmd_inv;
	gsl_spline *spline_bmd_inv;
};

class bmd_analysis
{
public:
	Eigen::MatrixXd MAP_ESTIMATE;
	Eigen::MatrixXd COV;
	bmd_cdf BMD_CDF;
	bool isExtra;
	double BMR;
	double MAP_BMD;
	double alpha;
	double MAP;
	contbmd type;
	std::vector<double> expected;

	bmd_analysis() : MAP_ESTIMATE(), COV(), BMD_CDF()
	{
	}
	bmd_analysis(const bmd_analysis &M)
	{
		BMD_CDF = M.BMD_CDF;
		MAP_ESTIMATE = M.MAP_ESTIMATE;
		MAP = M.MAP;
		COV = M.COV;
		BMR = M.BMR;
		MAP_BMD = M.MAP_BMD;
		alpha = M.alpha;
		type = M.type;
		expected = M.expected;
	}

	bmd_analysis &operator=(const bmd_analysis &M)
	{
		BMD_CDF = M.BMD_CDF;
		MAP_ESTIMATE = M.MAP_ESTIMATE;
		MAP = M.MAP;
		COV = M.COV;
		BMR = M.BMR;
		MAP_BMD = M.MAP_BMD;
		alpha = M.alpha;
		type = M.type;
		expected = M.expected;
		return *this;
	}
};

Eigen::MatrixXd convertresult_to_probs(Eigen::MatrixXd data);
double find_maBMDquantile(double q, std::vector<double> probs, std::list<bmd_analysis> analyses);
double ma_cdf(double dose, std::vector<double> probs, std::list<bmd_analysis> analyses);

/**********************************************************************
 * function bmd_analysis:
 * 		LL        -    likelihood previously defined
 *      PR        -    Prior
 * 	    fixedB    -    boolean vector of which parameters are fixed
 * 		fixedV    -    vector of what value parameters are fixed too.
 * 	    BMDType   -    The type of BMD we are searching for
 *      BMRF      -    Benchmark Response
 *      tail_prob -    The tail probability for the Hybrid approach
 *      isIncreasing - true if the BMD should be calculated for an increasing data set
 *      alpha     -    Defines the the lower and upper stopping point for integeration
 *      step_size -    The approximate step size
 *  returns:
 *      bmd_analysis class.
 * *******************************************************************/
template <class LL, class PR>
bmd_analysis bmd_analysis_CNC(LL likelihood, PR prior,
							  std::vector<bool> fixedB, std::vector<double> fixedV,
							  contbmd BMDType, double BMRF, double tail_prob,
							  bool isIncreasing, double alpha, double step_size,
							  Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1, 1))
{
	// value to return

	bmd_analysis rVal;
	// create the Continuous BMD models
	cBMDModel<LL, PR> model(likelihood, prior, fixedB, fixedV, isIncreasing);
	// Find the maximum a-posteriori and compute the BMD

	optimizationResult OptRes = findMAP<LL, PR>(&model, init);
	// DEBUG_OPEN_LOG("bmds.log", file);
	// DEBUG_LOG(file, "After findMap, optres= " << OptRes.result << ", MAP= " << OptRes.functionV << ", max_parms=\n" << OptRes.max_parms << "\n");

	double BMD = model.returnBMD(BMDType, BMRF, tail_prob);


	Eigen::MatrixXd result;
	std::vector<double> x;
	std::vector<double> y;

	if (!std::isinf(BMD) && !std::isnan(BMD))
	{
		int i = 0;
		do
		{
			result = profile_cBMDNC<LL, PR>(&model,
											BMDType,
											BMD, BMRF, tail_prob,
											step_size,
											(gsl_cdf_chisq_Pinv(1.0 - 2 * alpha, 1) + 0.1) / 2.0, // the plus 0.1 is to go slightly beyond
											isIncreasing);
			// if we have too few rows we decrease
			// the step size
			if (result.rows() <= 5)
				step_size *= 0.5;
			i++;
		} while (result.rows() <= 5 && !(i > 4));

		// Prepare the results to form an X/Y tuple
		// X - BMD values
		// Y - Cumulative Probabilities associated with the corresponding X row

		result = convertresult_to_probs(result);
		x.resize(result.rows());
		y.resize(result.rows());
	}


	// compute the CDF for the BMD posterior approximation
	if (!std::isinf(BMD) && !isnan(BMD) && BMD > 0 // flag numerical thins so it doesn't blow up.
		&& result.rows() > 5)
	{
		for (size_t i = 0; i < x.size(); i++)
		{
			x[i] = result(i, 0);
			y[i] = result(i, 1);
		}
		bmd_cdf cdf(x, y);
		rVal.BMD_CDF = cdf;
	}

	Eigen::MatrixXd estimated = model.log_likelihood.mean(OptRes.max_parms);

	rVal.expected.resize(estimated.rows());

	for (size_t i = 0; i < rVal.expected.size(); i++)
	{
		rVal.expected[i] = estimated(i, 0);
	}

	// std::cout << OptRes.max_parms << endl;
	// std::cout << model.log_likelihood.mean(OptRes.max_parms) << endl;
	// std::cout << model.log_likelihood.negLogLikelihood(OptRes.max_parms) << endl;

	rVal.MAP_BMD = BMD;
	rVal.BMR = BMRF;
	rVal.isExtra = false;
	rVal.type = BMDType;
	rVal.COV = model.varMatrix(OptRes.max_parms);
	rVal.MAP_ESTIMATE = OptRes.max_parms;
	// cout << "******&&&&& OptRes.functionV= " << OptRes.functionV << endl;
	rVal.MAP = OptRes.functionV;

	return rVal;
}

/**********************************************************************
 * function bmd_continuous_optimization:
 * 		LL        -    likelihood previously defined
 *      PR        -    Prior
 * 	    fixedB    -    boolean vector of which parameters are fixed
 * 	  	fixedV    -    vector of what value parameters are fixed too.
 *      isIncreasing - true if the BMD should be calculated for an increasing data set
 *  returns:
 *      Eigen::MatrixXd rVal <- the value that maximizes the model with the data.
 * *******************************************************************/
template <class LL, class PR>
Eigen::MatrixXd bmd_continuous_optimization(Eigen::MatrixXd Y, Eigen::MatrixXd X,
											Eigen::MatrixXd prior,
											std::vector<bool> fixedB, std::vector<double> fixedV,
											bool is_const_var,
											bool is_increasing,
											Eigen::MatrixXd init = Eigen::MatrixXd::Zero(10, 10))
{

	// value to return
	bool suff_stat = (Y.cols() == 3); // it is a SS model if there are three parameters

	LL likelihood(Y, X, suff_stat, is_const_var, is_increasing);
	PR model_prior(prior);
	Eigen::MatrixXd rVal;
	// create the Continuous BMD model

	cBMDModel<LL, PR> model(likelihood, model_prior, fixedB, fixedV, is_increasing);
	optimizationResult OptRes;
	// Find the maximum a-posteriori and compute the BMD

	if (init.cols() == 10 && init.rows() == 10)
	{
		OptRes = findMAP<LL, PR>(&model);
	}
	else
	{
		OptRes = findMAP<LL, PR>(&model, init);
	}
	rVal = OptRes.max_parms;
	return rVal;
}

/**********************************************************************
 * function bmd_continuous_optimization:
 * 		LL        -    likelihood previously defined
 *      PR        -    Prior
 * 	    fixedB    -    boolean vector of which parameters are fixed
 * 	  	fixedV    -    vector of what value parameters are fixed too.
 *      isIncreasing - true if the BMD should be calculated for an increasing data set
 *      degree    - the degree of the polynomial
 *  returns:
 *      Eigen::MatrixXd rVal <- the value that maximizes the model with the data.
 * *******************************************************************/
template <class LL, class PR>
Eigen::MatrixXd bmd_continuous_optimization(Eigen::MatrixXd Y, Eigen::MatrixXd X,
											Eigen::MatrixXd prior,
											std::vector<bool> fixedB, std::vector<double> fixedV,
											bool is_const_var,
											bool is_increasing,
											int degree)
{

	// value to return
	bool suff_stat = (Y.cols() == 3); // it is a SS model if there are three parameters

	LL likelihood(Y, X, suff_stat, is_const_var, degree); // for polynomial models
	PR model_prior(prior);
	Eigen::MatrixXd rVal;
	// create the Continuous BMD model
	cBMDModel<LL, PR> model(likelihood, model_prior, fixedB, fixedV, is_increasing);
	// Find the maximum a-posteriori and compute the BMD
	optimizationResult OptRes = findMAP<LL, PR>(&model);

	rVal = OptRes.max_parms;

	return rVal;
}

/**********
 * Struct for lambda function
 *
 */
template <class LL, class PR>
class fastBMDData
{
public:
	cBMDModel<LL, PR> *model;
	contbmd BMDType;
	double BMRF;
	double advP;
};
/**********************************************************
 * bmd_fast_BMD_cont:
 * This particular function computes BMD confidence intervals using Wald based
 * estimates vs. the standard profile likelihood.  It is primarily used in
 * high-throughput genomic DR analyses.
 ***********************************************************/
template <class LL, class PR>
bmd_analysis bmd_fast_BMD_cont(LL likelihood, PR prior,
							   std::vector<bool> fixedB, std::vector<double> fixedV,
							   contbmd BMDType, double BMRF, double tail_prob,
							   bool is_increasing,
							   Eigen::MatrixXd init = Eigen::MatrixXd::Zero(10, 10))
{

	bmd_analysis rVal;														   // Return Value
	optimizationResult oR;													   // Optimization result
	cBMDModel<LL, PR> model(likelihood, prior, fixedB, fixedV, is_increasing); // Model specification

	/************************************************
	 * Lambda Function
	 * for BMD Gradient computation
	 ************************************************/
	auto bmd = [](Eigen::MatrixXd a, void *b)
	{
		fastBMDData<LL, PR> *data = (fastBMDData<LL, PR> *)b;
		return data->model->returnBMD(a, data->BMDType, data->BMRF,
									  data->advP);
	};
	/*******************************************/

	if (init.rows() == 10 && init.cols() == 10) // Optimize
		oR = findMAP<LL, PR>(&model);
	else
		oR = findMAP<LL, PR>(&model, init, OPTIM_NO_FLAGS);

	Eigen::MatrixXd parms = oR.max_parms;

	/*
	 * Start Computing the BMD CI
	 */
	fastBMDData<LL, PR> data;
	data.model = &model;
	data.advP = tail_prob;
	data.BMDType = BMDType;
	data.BMRF = BMRF;

	double BMD = model.returnBMD(BMDType, BMRF, tail_prob);

	/*
	 * Calculate the Gradient for the delta Method
	 */

	double *g = new double[parms.rows()];

	gradient(parms, g, &data, bmd); // get the gradient vector
	Eigen::MatrixXd grad = parms * 0;
	for (int i = 0; i < grad.rows(); i++)
	{
		grad(i, 0) = g[i];
	}

	rVal.COV = model.varMatrix(parms);

	Eigen::MatrixXd var = grad.transpose() * rVal.COV * grad; // Delta Method Variance

	if (var(0) > 1e4)
	{
		var(0) = 1e4;
	} // If it is too large there are problems numerically
	  // make sure it isn't that large, but still large
	/*
	 * Compute CI using the Gaussian distribution. Here the Delta method is used again
	 * as the log(BMD) CI is computed.  This gaurantees the BMDL > 0.  It is then exponentiated.
	 */
	std::vector<double> x(500);
	std::vector<double> y(500);
	if (isnormal(var(0, 0)) && (var(0, 0) > 1e-7) && isnormal(log(BMD)))
	{

		for (size_t i = 0; i < x.size(); i++)
		{
			x[i] = double(i) / double(x.size());
			double q = x[i];
			y[i] = exp(gsl_cdf_gaussian_Pinv(q, (1.0 / BMD) * pow(var(0, 0), 0.5)) + log(BMD));
		}

		for (size_t i = y.size() - 1; i >= 1; i--)
		{
			if (y[i] == y[i - 1] || isinf(y[i]))
			{
				auto p1 = std::next(x.begin(), i);
				auto p2 = std::next(y.begin(), i);
				y.erase(p2);
				x.erase(p1);
				i = y.size() - 1;
			}
		}
	}
	else
	{
		x.resize(2);
		y.resize(2);
		x[0] = 0.0;
		x[1] = 1.0; // std::numeric_limits<double>::infinity();
		y[0] = 0.0;
		y[1] = 1.0;
	}

	if (isnormal(BMD) && BMD > 0.0 && // flag numerical thins so it doesn't blow up.
		x.size() > 6)
	{

		// fix dumb GSL numerical issues
		for (size_t i = 1; i < x.size(); i++)
		{
			if (x[i] <= x[i - 1])
			{
				for (size_t kk = i; kk < x.size(); kk++)
				{
					x[kk] = x[kk - 1] + 1e-6;
				}
			}
		}
		bmd_cdf cdf(x, y);
		rVal.BMD_CDF = cdf;
	}

	Eigen::MatrixXd estimated = model.log_likelihood.mean(oR.max_parms);

	rVal.expected.resize(estimated.rows());

	for (size_t i = 0; i < rVal.expected.size(); i++)
	{
		rVal.expected[i] = estimated(i, 0);
	}

	rVal.MAP_BMD = BMD;
	rVal.BMR = BMRF;
	rVal.isExtra = false;
	rVal.type = BMDType;
	rVal.MAP_ESTIMATE = oR.max_parms;
	rVal.MAP = oR.functionV;
	delete[] g;
	return rVal;
}

/**********************************************************************
 *  function: bmd_analysis_DNC(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
 *							 double BMR, bool isExtra, double alpha, double step_size)
 *
 * *******************************************************************/
template <class LL, class PR>
bmd_analysis bmd_analysis_DNC(Eigen::MatrixXd Y, Eigen::MatrixXd D, Eigen::MatrixXd prior,
							  std::vector<bool> fixedB, std::vector<double> fixedV, int degree,
							  double BMR, bool isExtra, double alpha, double step_size)
{

	LL dichotimousM(Y, D, degree);
	PR model_prior(prior);

	dBMDModel<LL, PR> model(dichotimousM, model_prior, fixedB, fixedV);
	signed int flags = OPTIM_USE_GENETIC | OPTIM_USE_SUBPLX;
	optimizationResult oR = findMAP<LL, PR>(&model, flags);

	bmd_analysis rVal;
	double BMD = isExtra ? model.extra_riskBMDNC(BMR) : model.added_riskBMDNC(BMR);

	Eigen::MatrixXd result;
	std::vector<double> x;
	std::vector<double> y;
	if (!std::isinf(BMD) && !std::isnan(BMD))
	{
		int i = 0;
		do
		{

			result = profile_BMDNC<LL, PR>(&model, isExtra,
										   BMD, BMR,
										   step_size, (gsl_cdf_chisq_Pinv(1.0 - 2 * alpha, 1) + 0.1) / 2.0, true);

			if (result.rows() <= 5)
				step_size *= 0.5;
			i++;
		} while (result.rows() <= 5 && !(i > 4));

		// Prepare the results to form an X/Y tuple
		// X - BMD values
		// Y - Cumulative Probabilities associated with the corresponding X row
		result = convertresult_to_probs(result);
		x.clear();
		y.clear();
		// fix the cdf so things don't implode
		for (int i = 0; i < result.rows(); i++)
		{
			if (!isnan(result(i, 0)) && !isinf(result(i, 0)))
			{
				y.push_back(result(i, 1));
				x.push_back(result(i, 0));
			}
		}
		// fix numerical quantile issues
		// so gsl doesn't go BOOM
		for (size_t i = 1; i < x.size(); i++)
		{
			if (x[i] <= x[i - 1])
			{
				for (size_t kk = i; kk < x.size(); kk++)
				{
					x[kk] = x[kk - 1] + 1e-6;
				}
			}
		}
	}

	if (!std::isinf(BMD) && !isnan(BMD) && BMD > 0 // flag numerical thins so it doesn't blow up.
		&& result.rows() > 5)
	{

		bmd_cdf cdf(x, y);
		rVal.BMD_CDF = cdf;
	}

	Eigen::MatrixXd estimated_p = model.log_likelihood.mean(oR.max_parms);
	rVal.expected.resize(estimated_p.rows());
	for (size_t i = 0; i < rVal.expected.size(); i++)
	{
		rVal.expected[i] = estimated_p(i, 0) * Y(i, 1);
	}

	rVal.MAP_BMD = BMD;
	rVal.BMR = BMR;
	rVal.isExtra = isExtra;
	rVal.COV = model.varMatrix(oR.max_parms);
	rVal.MAP_ESTIMATE = oR.max_parms;
	rVal.MAP = oR.functionV;

	return rVal;
}

template <class PR>
void RescaleContinuousModel(cont_model CM, Eigen::MatrixXd *prior, Eigen::MatrixXd *betas,
							double max_dose, double divisor,
							bool is_increasing, bool is_logNormal, bool is_const_var)
{

	Eigen::MatrixXd te_b = *betas;

	PR model_prior(*prior);
	divisor = max(1.0, divisor);

	int degree = 2;
	if (CM == cont_model::polynomial)
	{
		degree = prior->rows() - 2;
		if (!is_const_var)
		{
			degree--;
		}
	}
	Eigen::MatrixXd temp = rescale_parms(*betas, CM, max_dose, divisor, is_logNormal, degree);
	// fixme: in the future we might need to change a few things
	//  if there are more complicated priors
	int adverseR = 0;
	int nparms = te_b.rows();
	int tot_e = 1;

	switch (CM)
	{
	case cont_model::polynomial:

		if (!is_const_var)
		{
			tot_e = 2;
		}
		model_prior.scale_prior(divisor, 0);
		for (int i = 1; i < nparms - tot_e; i++)
		{
			model_prior.scale_prior(divisor, i);
			model_prior.scale_prior(pow(1 / max_dose, i), i);
		}
		break;
	case cont_model::funl:
		// b <- A[1] + A[2]*exp((doses-A[5])^2*(-A[6]))*(1/(1+exp(-(doses-A[3])/A[4])))

		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(divisor, 1);
		model_prior.scale_prior(max_dose, 2);
		model_prior.scale_prior(max_dose, 3);
		model_prior.scale_prior(max_dose, 4);
		//  model_prior.scale_prior(1/max_dose,5);
		// model_prior.add_mean_prior(1/max_dose*1/max_dose,5);
		if (!is_logNormal)
		{
			if (is_const_var)
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 5);
			}
			else
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 6);
			}
		}

		break;
	case cont_model::hill:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(divisor, 1);
		// model_prior.scale_prior(1/max_dose,1);
		model_prior.scale_prior(max_dose, 2);
		if (!is_logNormal)
		{
			if (is_const_var)
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 4);
			}
			else
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 5);
			}
		}

		break;
	case cont_model::hill_aerts:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(max_dose, 1);
		if (!is_logNormal){
			model_prior.add_mean_prior(2.0 * log(divisor), te_b.rows() - 1);
		}
		// else{
		// 	model_prior.scale_prior(log(divisor) / divisor, 0);
		// }
		break;
	case cont_model::lognormal_aerts: case cont_model::logskew_aerts: case cont_model::exp_aerts:
	case cont_model::gamma_aerts: case cont_model::invlomax_aerts:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(pow(max_dose, -te_b(3,0)), 1);
		if (!is_logNormal){
			model_prior.add_mean_prior(2.0 * log(divisor), te_b.rows() - 1);
		}// }else{
		// 	model_prior.scale_prior(log(divisor) / divisor, 0);
		// }
		break;
	case cont_model::gamma_efsa:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(pow(max_dose, -1.0), 1);
		if (!is_logNormal){
			model_prior.add_mean_prior(2.0 * log(divisor), te_b.rows() - 1);
		}
		// else{
		// 	model_prior.scale_prior(log(divisor) / divisor, 0);
		// }
		break;
	case cont_model::invlogskew_aerts: case cont_model::invgamma_aerts:
	case cont_model::invexp_aerts: case cont_model::lomax_aerts:// case cont_model::hill_aerts:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(pow(max_dose, te_b(3,0)), 1);
		if (!is_logNormal){
			model_prior.add_mean_prior(2.0 * log(divisor), te_b.rows() - 1);
		}
		// else{
		// 	model_prior.scale_prior(log(divisor) / divisor, 0);
		// }
		break;
	case cont_model::logistic_aerts: case cont_model::probit_aerts:
		model_prior.scale_prior(divisor, 2);
		model_prior.scale_prior(pow(max_dose, -te_b(3,0)), 1);
		if (!is_logNormal){
			model_prior.add_mean_prior(2.0 * log(divisor), te_b.rows() - 1);
		}
		// else{
		// 	model_prior.scale_prior(log(divisor) / divisor, 2);
		// }
		break;
	case cont_model::LMS:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(pow(max_dose, -1.0), 1);
		model_prior.scale_prior(pow(max_dose, -2.0), 3);
		if (!is_logNormal){
			model_prior.add_mean_prior(2.0 * log(divisor), te_b.rows() - 1);
		}
		// else{
		// 	model_prior.scale_prior(log(divisor) / divisor, 0);
		// }
		break;
	case cont_model::exp_3:
		adverseR = is_increasing ? NORMAL_EXP3_UP : NORMAL_EXP3_DOWN;
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(1 / max_dose, 1);
		if (!is_logNormal)
		{
			if (is_const_var)
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 4);
			}
			else
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 5);
			}
		}
		break;
	case cont_model::exp_5:
		adverseR = is_increasing ? NORMAL_EXP5_UP : NORMAL_EXP5_DOWN;
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(1 / max_dose, 1);
		if (!is_logNormal)
		{
			if (is_const_var)
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 4);
			}
			else
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 5);
			}
		}
		break;
	case cont_model::power:
		model_prior.scale_prior(divisor, 0);
		model_prior.scale_prior(divisor * pow(1 / max_dose, te_b(2, 0)), 1);
		// model_prior.scale_prior( divisor*(1/max_dose),2);

		if (!is_logNormal)
		{
			if (is_const_var)
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 3);
			}
			else
			{
				model_prior.add_mean_prior(2.0 * log(divisor), 4);
			}
		}
		break;

	default:
		break;
	}

	Eigen::MatrixXd temp2 = model_prior.get_prior();
	*prior = temp2;
	*betas = temp;
	return;
}

#endif
