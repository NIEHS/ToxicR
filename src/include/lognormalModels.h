// File: lognormalModels.h
// Purpose: 
// Implements the normal model likelihoods that are used for BMD analysis
//
//  Creator:       Matt Wheeler
//  Creation Date: 4/17/2018
#include "lognormal_likelihoods.h"
#include "cmodeldefs.h"

#ifndef  lognormalModels_H
#define lognormalModels_H

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "normalModels.h"
#include "cmodeldefs.h"

/*class loglognormalLLModel
 * This class defines a log-normal log-likelihood where the data 
 * basic assumption is log(y) ~ N(b0 + b1*d, sigma^2) 
 * 
 */

/*implements the basic model model Y = X*/
class lognormalLLModel : public lognormalLL {

public: 
	lognormalLLModel() {
	};


	lognormalLLModel(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
					          bool SS, bool CV = true) : lognormalLL(tY, tX,SS) {
		  
	};
	
	virtual cont_model mean_type(){
	  return cont_model::generic; 
	}	
	
	int    nParms() { 
				return 3; // linear regression + constant variance	
	}
	
	// BASIC MEAN FUNCTIONS
  virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta); 					    	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d);  

	//VARIANCE FUNCTIONS
	Eigen::MatrixXd variance(Eigen::MatrixXd theta); 	
	Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d); 


	//BASIC BMD Computation absolute to hybrid
	virtual double bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){return 0.0;}; 
  virtual double bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,double BPROB){return 0.0;};
    
    //BASIC BMD Computation absolute to hybrid
	virtual double bmd_absolute_bound(Eigen::MatrixXd theta,	 double BMD,double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_stdev_bound(Eigen::MatrixXd theta, 	     double BMD,double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_reldev_bound(Eigen::MatrixXd theta,       double BMD,double BMRF, bool isIncreasing){return 0.0;}; 
	virtual double bmd_point_bound(Eigen::MatrixXd theta,        double BMD,double BMRF, bool isIncreasing){return 0.0;};; 
	virtual double bmd_extra_bound(Eigen::MatrixXd theta,        double BMD,double BMRF, bool isIncreasing){return 0.0;} 
    virtual double bmd_hybrid_extra_bound(Eigen::MatrixXd theta, double BMD,double BMRF, bool isIncreasing,
										  double TAIL_PROB){return 0.0;};

	//Profiling requires good starting values for each of the individual models
	// these function do that for each type of risk
	virtual double bmd_start_extra_hybrid(unsigned n,
												   const double *b,
												   double *grad,
												   void   *data)  
	{
		return 0.0; 
	}
	////////////////////////////////////////////////////////
	virtual double bmd_start_absolute(unsigned n,
									  const double *b,
									  double *grad,
									  void   *data)  
	{
		return 0.0; 
	}

	virtual std::vector<double>  bmd_start_absolute_clean(std::vector<double> x,
																double BMRF, double BMD,
																bool isIncreasing)
	{
		return x;
	}
	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	virtual double bmd_start_stddev(unsigned n,
									const double *b,
									double *grad,
									void   *data) {
			return 0.0;
	};

	virtual std::vector<double> bmd_start_stddev_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing) {
			return x;
	};
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////

	virtual double bmd_start_reldev(unsigned n,
											 const double *b,
											 double *grad,
											 void   *data)  
	{
		return 0.0;
	}; 

	virtual std::vector<double>  bmd_start_reldev_clean(std::vector<double> x,
		double BMRF, double BMD,
		bool isIncreasing)
	{
		return x;
	};


	virtual double bmd_start_extra(unsigned n,
									const double *b,
									double *grad,
									void   *data) {
			return 0.0;
	};
	virtual std::vector<double> bmd_start_extra_clean(std::vector<double> x,
														double BMRF,
														double BMD,
														bool isIncreasing) {
		return x;
	};
    //////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	virtual double bmd_start_point(unsigned n,
									const double *b,
									double *grad,
									void   *data) {	
		return 0.0; 
	}

	virtual std::vector<double> bmd_start_point_clean(std::vector<double> x,
													double BMRF,
													double BMD,
													bool isIncreasing) {
		return x; 
	}
	
	//////////////////////////////////////////////////////
	virtual int type_of_profile(contbmd TYPE);
	virtual int parameter_to_remove(contbmd TYPE);
	//////////////////////////////////////////////////////
	virtual double  bmd_start_hybrid_extra(unsigned n,
													  const double *b,
													  double *grad,
													  void   *data){
		return 0.0; 
	}												  
	std::vector<double> bmd_start_hybrid_extra_clean(std::vector<double> x,
						  							 double BMRF, double BMD,
													 bool isIncreasing, 
													 double tail_prob){
		return x; 
	}	
	
	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
    Eigen::MatrixXd eqConst_gradient(Eigen::MatrixXd theta, contbmd type, 
									double BMD, double BMRF, bool isIncreasing, 
								    double tail_prob);
    ////////////////////	
	virtual double equality_boundG(Eigen::MatrixXd theta, contbmd BMDType,
									double BMD, double BMRF, bool isInc,
									double tail_prob);
	////////////////// Provide a starting value for the equality bound
	Eigen::MatrixXd starting_value(Eigen::MatrixXd theta, contbmd BMDType,
									double BMD, double BMRF, bool isInc,
							        double tail_prob,
							        std::vector<double> lb,
									std::vector<double> ub);





private: 


}; 

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

struct start_dataLN {
	lognormalLLModel *M; 
	Eigen::MatrixXd theta;
	double BMD;
	double BMRF;
	contbmd BMDType; 
	bool isIncreasing;
	double tail_prob;
};
#endif
