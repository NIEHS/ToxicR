#include <DichQlinearBMD_NC.h>

#ifdef R_COMPILATION
        //necessary things to run in R
        #include <RcppGSL.h>
        #include <RcppEigen.h>
#else
        #include <Eigen/Dense>
#endif


double QLINEAR_BMD_EXTRA_NC_INEQUALITY(Eigen::MatrixXd theta, void* data) {
	log_qlinear_inequality *M = (log_qlinear_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double g = QLINEAR_G(theta(0, 0));
	double Z = QLINEAR_EXTRA_Z(g, BMR); //note BMD is a placeholder
	Z = Z / BMD;
	double rV = 0.0;
	rV = (geq) ? inequality - Z : Z - inequality;
	return rV;
}

double QLINEAR_BMD_ADDED_NC_INEQUALITY(Eigen::MatrixXd theta, void* data) {
	log_qlinear_inequality *M = (log_qlinear_inequality*)data;
	double inequality = M->inequality;
	double BMD = M->BMD;
	double BMR = M->BMR;
	bool   geq = M->geq;

	double g = QLINEAR_G(theta(0, 0));

	double Z = QLINEAR_ADDED_Z(g, BMR);
	Z = Z/BMD;
	double rV = 0.0;

	rV = (geq) ? inequality - Z : Z - inequality;

	return rV;
}
