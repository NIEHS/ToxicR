#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <cfloat>
#include <nlopt.hpp>

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
  
#endif

#include "bmd_calculate.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std; 

Eigen::MatrixXd convertresult_to_probs(Eigen::MatrixXd data) {
	Eigen::MatrixXd result = data;
	bool is_afterMax = false;
	for (int i = 0; i < result.rows(); i++) {
		if (result(i, 2) == 666) {
			is_afterMax = true;
		}
		result(i,0 ) = result(i, 0) < 0 ?  1e4: result(i, 0); // numerical zero for calculations
		result(i, 0) = is_afterMax ? .5 + gsl_cdf_chisq_P(2.0*result(i, 0), 1.0) / 2.0 : .5 - gsl_cdf_chisq_P(2.0*result(i, 0), 1.0) / 2.0;

		if (i > 0)
			result(i, 2) = result(i, 0) - result(i - 1,0);

	}
	// make sure the difference has some small probability of change. 
	for (int i = 1; i < result.rows(); i++) {
		if (result(i, 2) <= 0.0) {
			result(i, 0) = result(i - 1, 0) + 1e-4;
			if (i != result.rows() - 1) {                          // if it is not the last one we have 
				result(i+1, 2) = result(i+1, 0) - result(i,0);    // we must update the difference
			}
		}
	}

	// now we look at the individual rows and make sure there is no nan 
	// if it is a nan we remove that row. 
	int numBad = 0;

	for (int i = 0; i < result.rows(); i++) {
		if (isnan(result(i, 0)) || isinf(result(i, 0))) {
			result(i, 2) = 77;
			numBad++;
		}
	}

	if (numBad) {
		Eigen::MatrixXd tempMat(result.rows() - numBad, result.cols());
		int counter = 0;
		for (int i = 0; i < result.rows(); i++) {
			if (result(i, 2) != 77) {
				tempMat.row(counter) = result.row(i);
				tempMat(counter, 2) = data(i, 2);
				counter++;
			}
		}
		result = tempMat;
	}
	return result;
}




double ma_cdf(double dose, std::vector<double> probs, std::list<bmd_analysis> analyses) {
	double rVal = 0.0;
	double temp;
	int i = 0;
	for (bmd_analysis b : analyses) {
		temp = b.BMD_CDF.P(dose);
		rVal += probs[i] * temp; i++;
	}
	return rVal;
};

double find_maBMDquantile(double q, std::vector<double> probs, std::list<bmd_analysis> analyses) {
	double temp = 0.0;
	double max = 1.0; double min = 0.0;
	bool quantile_isNotBounded = true;
	while (quantile_isNotBounded && !std::isinf(temp) && !std::isnan(temp)) {
		temp = ma_cdf(max, probs, analyses);
		if (temp > q) {
			quantile_isNotBounded = false;
		}
		max *= 2;
	}

	double test = (max + min) / 2.0;
	temp = ma_cdf(test, probs, analyses);
	while ((fabs(log(temp / q)) > 1e-8) && !std::isinf(temp) && !std::isnan(temp)) { // use the log ratio in the calculation as it is invariant to the size of q
		if (temp > q) {
			max = test;
		}
		else {
			min = test;
		}
		test = (max + min) / 2.0;
		temp = ma_cdf(test, probs, analyses);
	}
    if (std::isinf(temp) || std::isnan(temp))
      return NAN;
    else
      return test;
}

