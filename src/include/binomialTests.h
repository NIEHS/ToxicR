#pragma once

#ifndef BINOMIAL_TESTS_H
#define BINOMIAL_TESTS_H

#ifdef R_COMPILATION
    //necessary things to run in R    
      #include <R.h>
    #include <Rmath.h>    
    #include <Rinternals.h>
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include "log_likelihoods.h"

/*implements the basic model model Y = X*/
class binomialLLTESTA1 : public binomialLL {

			public:
				binomialLLTESTA1(Eigen::MatrixXd tY, Eigen::MatrixXd tX) : binomialLL(tY, tX) {
					std::vector<double> vec(tX.data(), tX.data() + tX.rows() * tX.cols());
					std::sort(vec.begin(), vec.end());
					vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
					udoses = vec; // this should be the unique dose group

					meanX = Eigen::MatrixXd::Zero(tY.rows(), udoses.size());

					for (int i = 0; i < meanX.rows(); i++)
					{
						for (int j = 0; j < udoses.size(); j++) {
							meanX(i, j) = udoses[j] == X(i, 0) ? 1.0 : 0.0;
						}

					}


			};

				
			int    nParms() { return meanX.cols(); }; // the model is fully saturated
														  // if it is a sufficient statistics model
				
			virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
				return   mean(theta, meanX);
			};

			virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
				return variance(theta, meanX);
			};

			virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
				Eigen::MatrixXd p = 1 / (1 + exp(-(meanX*theta).array()));
				return  p.array()*(1 - p.array());
			};

			virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
				return   1 / (1 + exp(-(meanX*theta).array()));
			};

		// BASIC MEAN FUNCTIONS
		//		virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta);
		//		virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d);

				//VARIANCE FUNCTIONS
		//		Eigen::MatrixXd variance(Eigen::MatrixXd theta);
		//		Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d);

			public:
					std::vector<double> udoses; // unique dose groups 
					Eigen::MatrixXd     meanX; 

};

/*implements the basic model model Y = X*/
class binomialLLTESTA2 : public binomialLL {

public:
	binomialLLTESTA2(Eigen::MatrixXd tY, Eigen::MatrixXd tX) : binomialLL(tY, tX) {
	};
	
	int    nParms() { return 1; };
								  
	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return   mean(theta, X);
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return variance(theta, X);
	};

	virtual Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
		Eigen::MatrixXd ones = X;
		ones.setOnes();
		return   1 / (1 + exp(-(ones*theta(0, 0)).array()));
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd X) {
		Eigen::MatrixXd ones = X; 
		ones.setOnes(); 
		return   1 / (1 + exp(-(ones*theta(0,0)).array()));
	};

	// BASIC MEAN FUNCTIONS
	//		virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta);
	//		virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d);

	//VARIANCE FUNCTIONS
	//		Eigen::MatrixXd variance(Eigen::MatrixXd theta);
	//		Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d);

public:
	std::vector<double> udoses; // unique dose groups 
	Eigen::MatrixXd     meanX;

};

#endif
