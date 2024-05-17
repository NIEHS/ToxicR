#define STRICT_R_HEADERS

#include "DichLogProbitBMD_NC.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#ifdef R_COMPILATION
// necessary things to run in R
#ifdef ToxicR_DEBUG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <RcppEigen.h>
#pragma GCC diagnostic pop
#else
#include <RcppEigen.h>
#endif
#include <RcppGSL.h>
#else
#include <Eigen/Dense>
#endif

double logProbit_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void *data) {
  log_probit_inequality *M = (log_probit_inequality *)data;
  double inequality = M->inequality;
  double BMD = M->BMD;
  double BMR = M->BMR;
  bool geq = M->geq;

  double g = LOGPROBIT_G(theta(0, 0));
  double a = LOGPROBIT_A(theta(1, 0));
  double Z = LOGPROBIT_EXTRA_Z(g, a, BMR); // note BMD is a placeholder
  Z = Z / log(BMD);
  double rV = 0.0;
  rV = (geq) ? inequality - Z : Z - inequality;
  return rV;
}

double logProbit_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void *data) {
  log_probit_inequality *M = (log_probit_inequality *)data;
  double inequality = M->inequality;
  double BMD = M->BMD;
  double BMR = M->BMR;
  bool geq = M->geq;

  double g = LOGPROBIT_G(theta(0, 0));
  double a = LOGPROBIT_A(theta(1, 0));

  double Z = LOGPROBIT_ADDED_Z(g, a, BMR);
  Z = Z / log(BMD);
  double rV = 0.0;

  rV = (geq) ? inequality - Z : Z - inequality;

  return rV;
}
