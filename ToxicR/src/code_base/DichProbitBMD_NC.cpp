#include <DichProbitBMD_NC.h>

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

double PROBIT_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data) {
	probit_inequality *M = (probit_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double a = PROBIT_A(theta(0, 0));
	double Z = PROBIT_EXTRA_Z(a, BMR); //note BMD is a placeholder
	Z = Z / BMD;
	double rV = 0.0;
	rV = (geq) ? inequality - Z : Z - inequality;
	return rV;
}

double PROBIT_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data) {
    probit_inequality *M = (probit_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double a = PROBIT_A(theta(0, 0));
	double Z = PROBIT_ADDED_Z(a, BMR);
	Z = pow(Z, a) / pow(BMD, a);
	double rV = 0.0;

	rV = (geq) ? inequality - Z : Z - inequality;

	return rV;
}
