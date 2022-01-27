#pragma once
#ifndef lognormal_POLYNOMIAL_NCH
#define lognormal_POLYNOMIAL_NCH


#ifdef R_COMPILATION
    //necessary things to run in R
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else
    #include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>
#include "log_likelihoods.h"
#include "lognormalModels.h"

/* class lognormalPOLYNOMIAL_BMD_NC : public LL
 *
 * This class defines a lognormal log-likelihood where
 * the POLYNOMIAL model is the model that is fit to the data.
 * NC stands for no covariates, for future extentions that
 * may include covariates.
 */

/*implements the basic model model Y = X*/
class lognormalPOLYNOMIAL_BMD_NC : public lognormalLLModel {
	public:

	lognormalPOLYNOMIAL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS,
								 int degree) : lognormalLLModel(tY, tX,SS) {
			deg = degree;
	};

	lognormalPOLYNOMIAL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS, bool CV,
                              int degree) : lognormalLLModel(tY, tX,SS) {
	    deg = degree;
	};
	lognormalPOLYNOMIAL_BMD_NC(){

	};

	virtual cont_model mean_type(){
	  return cont_model::polynomial; 
	}	
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	}

	virtual int type_of_profile(contbmd TYPE);
	virtual int parameter_to_remove(contbmd TYPE);

	lognormalPOLYNOMIAL_BMD_NC(lognormalPOLYNOMIAL_BMD_NC &M){
		 sufficient_statistics = M.sufficient_statistics;
		 Y = M.Y;
		 X = M.X;
		 deg = M.deg;
	};

	int    nParms() {
				return 2+deg; // Power model regression + constant variance

	}




	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta,Eigen::MatrixXd d){

		Eigen::MatrixXd rV = theta(0, 0) + 0.0*pow(d.array(), 0);
		Eigen::MatrixXd temp;
		for (int i = 1; i < deg + 1; i++) {
			 temp  =   theta(i, 0)*pow(d.array(), double(i)); // sum up each degree of the polynomial
			 rV.col(0) += temp;
		}
		Eigen::MatrixXd returnV = log(rV.array());  // log the median

		return returnV;

	}

	// return true if it is a increasing function

	//BASIC BMD Computation absolute to hybrid
	double bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing);
	double bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing);
	double bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing);
	double bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing);
	double bmd_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing);
  double bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,double BPROB);

    //BASIC BMD Computation absolute to hybrid
	double bmd_absolute_bound(Eigen::MatrixXd theta,	 double BMD,double BMRF, bool isIncreasing);
	double bmd_stdev_bound(Eigen::MatrixXd theta, 	     double BMD,double BMRF, bool isIncreasing);
	double bmd_reldev_bound(Eigen::MatrixXd theta,       double BMD,double BMRF, bool isIncreasing);
	double bmd_point_bound(Eigen::MatrixXd theta,        double BMD,double BMRF, bool isIncreasing);
	double bmd_extra_bound(Eigen::MatrixXd theta,        double BMD,double BMRF, bool isIncreasing);
  double bmd_hybrid_extra_bound(Eigen::MatrixXd theta, double BMD,double BMRF, bool isIncreasing,
								  double TAIL_PROB);

	///
	/*Eigen::MatrixXd bmd_start_extra_hybrid(Eigen::MatrixXd theta, double BMD,
											double BMRF, bool isIncreasing,
											double TAIL_PROB);*/
	//
	virtual double bmd_start_absolute(unsigned n,
									  const double *b,
									  double *grad,
									  void   *data);
	virtual std::vector<double>  bmd_start_absolute_clean(std::vector<double> x,
															double BMRF, double BMD, bool isIncreasing);

	//
	/*Eigen::MatrixXd bmd_start_stdev(Eigen::MatrixXd theta, double BMD,
									double BMRF, bool isIncreasing,
									double TAIL_PROB);
	//*/
	virtual double bmd_start_reldev(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_reldev_clean(std::vector<double> x,
										  double BMRF,
										  double BMD,
										  bool isIncreasing);

	/////////////////////////////////////////////////////////////////////////
	// standard deviation starting value code
	//
	/////////////////////////////////////////////////////////////////////////
	virtual double bmd_start_stddev(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_stddev_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing);

	/////////////////////////////////////////////////////////////////////////
	// standard deviation starting value code
	/////////////////////////////////////////////////////////////////////////

	virtual double bmd_start_extra(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_extra_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing);

	/////////////////////////////////////////////////////////
	virtual double bmd_start_point(unsigned n,
									const double *b,
									double *grad,
									void   *data);

	virtual std::vector<double> bmd_start_point_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing);

	/////////////////////////////////////////////////////////

	virtual double  bmd_start_hybrid_extra(unsigned n,
													  const double *b,
													  double *grad,
													  void   *data);

	std::vector<double> bmd_start_hybrid_extra_clean(std::vector<double> x,
						  							 double BMRF, double BMD,
													 bool isIncreasing, double tail_prob);

	/////////////////////////////////////////////////////////
	/*
	Eigen::MatrixXd bmd_start_point(Eigen::MatrixXd theta, double BMD,
									double BMRF, bool isIncreasing,
									double TAIL_PROB);
	//
	Eigen::MatrixXd bmd_start_extra(Eigen::MatrixXd theta, double BMD,
									double BMRF, bool isIncreasing,
									double TAIL_PROB);*/

	private:
		int deg;
	// add BMD specific stuff here
};

#endif
