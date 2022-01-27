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

void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd,void*)> math_func) {

	Eigen::VectorXd h(v.rows());
	double mpres = pow(1.0e-16, 0.333333);
	double x, temp;
	double derEst;
	double f1 = 0.0; double f2 = 0.0;
	Eigen::MatrixXd  hvector = v;


	for (int i = 0; i < v.rows(); i++)
	{
		x = v(i, 0);
		if (fabs(x) > DBL_EPSILON) {
			h[i] = mpres * (fabs(x));
			temp = x + h[i];
			h[i] = temp - x;
		}
		else {
			h[i] = mpres;
		}

	}
	/*find the gradients for each of the variables in the likelihood*/
	for (int i = 0; i < v.rows(); i++) {
		/*perform a finite difference calculation on the specific derivative*/
		x = v(i, 0);
		// add h
		hvector(i, 0) = x + h[i];
		f1 = math_func(hvector,data);
		// subtract h
		hvector(i, 0) = x - h[i];
		f2 = math_func(hvector,data);
	
		// estimate the derivative
		g[i]  = (f1 - f2) / (2.0*h[i]);

		hvector(i, 0) = x;
		
	}
	return;
}

