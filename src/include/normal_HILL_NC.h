#pragma once

#ifndef normal_HILL_NCH
#define normal_HILL_NCH
#include "normalModels.h"
#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
    
#endif

#include <gsl/gsl_randist.h>


/*class normalHILL_BMD_NC : public LL
 * This class defines a normal log-likelihood where 
 * the HILL model is the model that is fit to the data.
 * NC stands for no covariates, for future extentions that 
 */

/*implements the basic model model Y = X*/



class normalHILL_BMD_NC : public normalLLModel {
	public: 

	normalHILL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS,
					          bool CV, int junk) : normalLLModel(tY, tX,SS,CV) {
		  // if it is a sufficient statistics model 
	};
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	}
	normalHILL_BMD_NC(){
	
	};
	
	virtual cont_model mean_type(){
	  return cont_model::hill; 
	}	
	
	normalHILL_BMD_NC(normalHILL_BMD_NC &M){
		 sufficient_statistics = M.sufficient_statistics; 
		 Y = M.Y; 
		 X = M.X; 
		 constant_variance = M.constant_variance; 
	};
	
	int    nParms() { 
				if (isConstVar()){
					return 5; // Hill regression + constant variance
				} else{
					return 6; // Hill regression + variance proportional to mean
				}  ; 		  
	}
	
	virtual int parameter_to_remove(contbmd TYPE); 

	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta,Eigen::MatrixXd d){
		double gamma = theta(0,0); 
		double nu = theta(1,0); 
		double k  = theta(2,0); 
		double n_exp = theta(3,0); 
		
		Eigen::MatrixXd rV = gamma + nu*pow(d.array(),n_exp)/(pow(k,n_exp)+pow(d.array(),n_exp));		
		
		return rV; 
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

//	virtual int type_of_profile(contbmd TYPE);


	private:
	// add BMD specific stuff here
}; 

#endif
