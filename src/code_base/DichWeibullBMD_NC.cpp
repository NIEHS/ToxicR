#include <DichWeibullBMD_NC.h>

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif

double WEIBULL_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data){
	log_weibull_inequality *M = (log_weibull_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double g = WEIBULL_G(theta(0, 0));
	double a = WEIBULL_A(theta(1, 0));
	double Z = WEIBULL_EXTRA_Z(g, a, BMR); //note BMD is a placeholder
	Z = pow(Z, a)/pow(BMD,a);
	double rV = 0.0;
	rV = (geq) ? inequality - Z : Z - inequality;

	return rV;
}

double WEIBULL_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data) {
	log_weibull_inequality *M = (log_weibull_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double g = WEIBULL_G(theta(0, 0));
	double a = WEIBULL_A(theta(1, 0));

	double Z = WEIBULL_ADDED_Z(g, a, BMR);
	Z = pow(Z,a) / pow(BMD,a);
	double rV = 0.0;

	rV = (geq)? inequality - Z: Z - inequality;

	return rV;
}

