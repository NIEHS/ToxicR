#include <DichGammaBMD_NC.h>

#ifdef R_COMPILATION
        //necessary things to run in R
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

double GAMMA_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data){
	log_gamma_inequality *M = (log_gamma_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double g = GAMMA_G(theta(0, 0));
	double a = GAMMA_A(theta(1, 0));
	double Z = GAMMA_EXTRA_Z(g, a, BMR); //note BMD is a placeholder
	Z = Z/BMD;
	double rV = 0.0;
	rV = (geq) ? inequality - Z : Z - inequality;
	return rV;
}

double GAMMA_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data) {
	log_gamma_inequality *M = (log_gamma_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double g = GAMMA_G(theta(0, 0));
	double a = GAMMA_A(theta(1, 0));

	double Z = GAMMA_ADDED_Z(g, a, BMR);
	Z = Z / BMD;
	double rV = 0.0;

	rV = (geq)? inequality - Z: Z - inequality;

	return rV;
}
