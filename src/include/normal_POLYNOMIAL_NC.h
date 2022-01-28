/*
 * Copyright 2020  NIEHS <matt.wheeler@nih.gov>
 * 
 *
 *Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
 *and associated documentation files (the "Software"), to deal in the Software without restriction, 
 *including without limitation the rights to use, copy, modify, merge, publish, distribute, 
 *sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
 *is furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all copies 
 *or substantial portions of the Software.
 
 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 *INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
 *PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
 *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
 *CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 *OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * 
 */
#pragma once
#ifndef normal_POLYNOMIAL_NCH
#define normal_POLYNOMIAL_NCH
#include "normalModels.h"
#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
    #include <gsl/gsl_randist.h>
#endif

/* class normalPOLYNOMIAL_BMD_NC : public LL
 *
 * This class defines a normal log-likelihood where 
 * the POLYNOMIAL model is the model that is fit to the data.
 * NC stands for no covariates, for future extentions that 
 * may include covariates.
 */

/*implements the basic model model Y = X*/
class normalPOLYNOMIAL_BMD_NC : public normalLLModel {
	public: 

	normalPOLYNOMIAL_BMD_NC(Eigen::MatrixXd tY, Eigen::MatrixXd tX,bool SS,
						 bool CV, int degree) : normalLLModel(tY, tX,SS,CV) {
			// if it is a sufficient statistics model 
			deg = degree; 
	};
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	}

	normalPOLYNOMIAL_BMD_NC(){
	
	};
	
	normalPOLYNOMIAL_BMD_NC(normalPOLYNOMIAL_BMD_NC &M){
		 sufficient_statistics = M.sufficient_statistics; 
		 Y = M.Y; 
		 X = M.X; 
		 deg = M.deg; 
		 constant_variance = M.constant_variance; 
	};
	
	int    nParms() { 
				if (isConstVar()){
					return 2+deg; // Power model regression + constant variance
				} else{
					return 3+deg; // Hill regression + variance proportional to mean
				}; 		  
	}
	
	virtual int type_of_profile(contbmd TYPE);
	virtual int parameter_to_remove(contbmd TYPE);

	virtual cont_model mean_type(){
	  return cont_model::polynomial; 
	}	
	
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta,Eigen::MatrixXd d){
				
		Eigen::MatrixXd rV = theta(0, 0) + 0.0*pow(d.array(), 0); 
		Eigen::MatrixXd temp; 
		
		for (int i = 1; i < deg + 1; i++) {
			 temp  =   theta(i, 0)*pow(d.array(), double(i)); // sum up each degree of the polynomial
			 rV.col(0) += temp; 
		}
		
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
