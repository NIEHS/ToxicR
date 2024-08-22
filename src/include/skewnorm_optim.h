#pragma once

#ifndef skewnorm_optimH
#define skewnorm_optimH
#include "owenst_asa076.h"
#ifdef R_COMPILATION
// necessary things to run in R
// necessary things to run in R
// #ifdef ToxicR_DEBUG
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wignored-attributes"
// #include <RcppEigen.h>
// #pragma GCC diagnostic pop
// #else
#include <RcppEigen.h>
// #endif
#include <RcppGSL.h>
#else
#include <Eigen/Dense>

#endif

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

struct skewnorm_params {
  double xi, q;
};

double skewnorm_pdf(double x, void *params);
double skewnorm_cdf(double x, void *params);
void skewnorm_fdf(double x, void *params, double *y, double *dy);
double skewnorm_quantile(double q, double xi);

#endif
