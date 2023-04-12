#include "stdafx.h" // Precompiled header - does nothing if building R version
#include "skewnorm_optim.h"
#include "owenst_asa076.h"

#ifdef R_COMPILATION
    //necessary things to run in R
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else
    #include <Eigen/Dense>

#endif

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <iostream>
#include <cmath>
#include <math.h>
using namespace std;

double skewnorm_pdf(double x, void *params){
	  struct skewnorm_params *p
	    = (struct skewnorm_params *) params;
	  double xi = p->xi;

	return 2.0*gsl_ran_ugaussian_pdf(x) * gsl_cdf_ugaussian_P(xi*x);
}

double skewnorm_cdfq(double x, void *params){
	  struct skewnorm_params *p
	    = (struct skewnorm_params *) params;
	  double xi = p->xi;
	  double q = p->q;
	return gsl_cdf_ugaussian_P(x) - 2.0*tfn(x, xi) - q;
}

void skewnorm_fdf(double x, void *params, double *y, double *dy){
	  struct skewnorm_params *p
	    = (struct skewnorm_params *) params;
	*y = skewnorm_cdfq(x, p);
	*dy = skewnorm_pdf(x, p);
}

double skewnorm_quantile(double q, double xi){

//	int status;
//	int iter = 0, max_iter = 5000;
//	const gsl_root_fdfsolver_type *T;
//	gsl_root_fdfsolver *s;
//	double x0 = 1.0, x = 0.0;
//	gsl_function_fdf FDF;
//
//	struct skewnorm_params params = {xi, q};
//
//
//	FDF.f = &skewnorm_cdfq;
//	FDF.df = &skewnorm_pdf;
//	FDF.fdf = &skewnorm_fdf;
//	FDF.params = &params;
//
//	T = gsl_root_fdfsolver_newton; // or steffenson might be faster (fdfsolver_steffenson)
//	s = gsl_root_fdfsolver_alloc (T);
//	gsl_root_fdfsolver_set (s, &FDF, x);
//
//	do
//	{
//		iter++;
//		status = gsl_root_fdfsolver_iterate (s);
//		x0 = x;
//		x = gsl_root_fdfsolver_root (s);
//		status = gsl_root_test_delta (x, x0, 0, 1e-3);
//	}
//	while (status == GSL_CONTINUE && iter < max_iter);
//
//	gsl_root_fdfsolver_free (s);
//	return x;

	  int status;
	  int iter = 0, max_iter = 1000;
	  const gsl_root_fsolver_type *T;
	  gsl_root_fsolver *s;
	  double r = 0;
	  double x_lo = -1e4, x_hi = 1e4;
	  gsl_function F;
	  struct skewnorm_params params = {xi, q};

	  F.function = &skewnorm_cdfq;
	  F.params = &params;

	  T = gsl_root_fsolver_brent;
	  s = gsl_root_fsolver_alloc (T);
	  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (s);
	      r = gsl_root_fsolver_root (s);
	      x_lo = gsl_root_fsolver_x_lower (s);
	      x_hi = gsl_root_fsolver_x_upper (s);
	      status = gsl_root_test_interval (x_lo, x_hi,
	                                       0, 0.0001);

	    }
	  while (status == GSL_CONTINUE && iter < max_iter);

	  gsl_root_fsolver_free (s);

	  return r;
}
