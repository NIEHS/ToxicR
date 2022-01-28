#pragma once

#ifndef normal_SKEW_NORMAL_NCH
#define normal_SKEW_NORMAL_NCH
#include "normalModels.h"
#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
    
#endif

#include <gsl/gsl_randist.h>


/*class normalFUNL_BMD_NC : public LL
 * This class defines a normal log-likelihood where 
 * the HILL model is the model that is fit to the data.  PROFILE_EQUALITY
 * NC stands for no covariates, for future extentions that 
 */

/*implements the basic model model Y = X*/
class normalFUNL_BMD_NC : public normalLLModel {
	public: 

	normalFUNL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS,
					          bool CV, int junk) : normalLLModel(tY, tX,SS,CV) {
		  // if it is a sufficient statistics model 
	};
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	}
	normalFUNL_BMD_NC(){
	
	};
	
	virtual cont_model mean_type(){
	  return cont_model::funl; 
	}	
	
	normalFUNL_BMD_NC(normalFUNL_BMD_NC &M){
		 sufficient_statistics = M.sufficient_statistics; 
		 Y = M.Y; 
		 X = M.X; 
		 constant_variance = M.constant_variance; 
	};
	
	int    nParms() { 
				if (isConstVar()){
					return 7; // Hill regression + constant variance
				} else{
					return 8; // Hill regression + variance proportional to mean
				}   		  
	}
	
	virtual int parameter_to_remove(contbmd TYPE); 
	virtual int type_of_profile(contbmd TYPE);
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta,Eigen::MatrixXd d){
		double b1 = theta(0,0); 
		double b2 = theta(1,0); 
		double b3 = theta(2,0); 
		double b4 = theta(3,0);
		double b5 = theta(4,0); 
		double b6 = theta(5,0); 
		
		Eigen::MatrixXd cdf = d; 
		Eigen::MatrixXd pdf = d; 
		for(int i = 0; i < d.rows(); i++){
		  cdf(i,0) =  1/(1+exp(-(d(i,0) - b3)/b4)); 
		  pdf(i,0) = exp(-exp(b6)*(d(i,0) - b5)*(d(i,0) - b5)); 
		}
		Eigen::MatrixXd rV = b1 + b2*pdf.array()*cdf.array();		
		
		return rV; 
	}


  double FUNL_mean(double d, Eigen::MatrixXd A){
    
    double rV = A(0,0) + A(1,0)*exp(-exp(A(5,0))*(d-A(4,0))*(d-A(4,0)))/(1  + exp(-(1/A(3,0))*(d-A(2,0))));
    return rV; 
  }
	
  double dFUNL_mean(double d,Eigen::MatrixXd A){
    double rV = -2*exp(A(5,0))*(d-A(4,0)) + (1/A(3,0))*exp(-(1/A(3,0))*(d-A(2,0)))/(1  + exp(-(1/A(3,0))*(d-A(2,0))));
    return rV; 
  }
	
  double d2FUNL_mean(double d,Eigen::MatrixXd A){
    double temp = (1/A(3,0)); temp *= temp; 
    double temp2= ((1  + exp(-(1/A(3,0))*(d-A(2,0))))); temp2 *= temp2; 
    
    double rV = -2*exp(A(5,0)) - (temp)*exp(-(1/A(3,0))*(d-A(2,0)))/temp2; 
    return rV; 
  }
	
	double findOptim(Eigen::MatrixXd A,bool isIncreasing){
	  double cv = 0.5; 
	  double nv = (cv - dFUNL_mean(cv,A)/d2FUNL_mean(cv,A))*.7; 
	  double test = nv - cv; 
	  int i = 0; 
	  
	  while (i < 250 && fabs(test) > 1e-8){
	    nv = cv - dFUNL_mean(cv,A)/d2FUNL_mean(cv,A); 
	    test = nv - cv; 
	    cv = nv; 
	    i++; 
	  }
	  
	  return cv; 
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
