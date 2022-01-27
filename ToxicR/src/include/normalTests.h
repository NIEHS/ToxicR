#pragma once
#ifndef normalTestsH
#define normalTestsH

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <R.h>
    #include <Rmath.h>    
    #include <Rinternals.h>
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
    #include <gsl/gsl_randist.h>
#endif

#include "cmodeldefs.h"

/*implements the basic model model Y = X*/
class normalLLTESTA1 : public normalLL {

			public:
				normalLLTESTA1() {};
				normalLLTESTA1(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
												   bool SS) : normalLL(tY, tX, SS){
					// if it is a sufficient statistics model
					std::vector<double> vec(tX.data(), tX.data() + tX.rows() * tX.cols());
					std::sort(vec.begin(), vec.end());
					vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
					udoses = vec; // this should be the unique dose group
					
					meanX = Eigen::MatrixXd::Zero(tY.rows(), udoses.size()); 
					
					for (int i = 0; i < meanX.rows(); i++)
					{
						for (int j = 0; j < udoses.size(); j++) {
							meanX(i, j) = udoses[j] == X(i, 0) ? 1.0 : 0.0; 
						}
						
					}
					

				};

				int    nParms() {
					return meanX.cols() + 1; 
				};

				virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
					return mean(theta, X);
				}; 

				virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
					Eigen::MatrixXd temp(theta.size() - 1, 1); 
					for (int i = 0; i < temp.size(); i++) {
						temp(i, 0) = theta(i, 0); 
					}
					Eigen::MatrixXd rV = meanX * temp; 
					return rV;
				}; 

				Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
					return variance(theta, X);
				}; 

				 Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
					double var = exp(theta(theta.rows() - 1, 0));
					Eigen::MatrixXd rV = Eigen::MatrixXd::Ones(d.rows(), 1)*var;
					return rV;
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

////////////////////////////////////////////////
class normalLLTESTA2 : public normalLL {

public:
	normalLLTESTA2() {};
	normalLLTESTA2(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
		bool SS) : normalLL(tY, tX, SS) {
		// if it is a sufficient statistics model
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

	// mean and variance for every dose group
	int    nParms() {
		return 2 * meanX.cols();
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		Eigen::MatrixXd temp(theta.size()/2, 1);
		for (int i = 0; i < temp.size(); i++) {
			temp(i, 0) = theta(i, 0);
		}
	
		Eigen::MatrixXd rV = meanX * temp;
		return rV;
	};

	Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return variance(theta, X);
	};

	Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		Eigen::MatrixXd temp(theta.size() / 2 , 1);
		for (int i = 0; i < temp.size(); i++) {
			temp(i, 0) = exp(theta(theta.size() / 2 + i, 0));
		}
		Eigen::MatrixXd rV = meanX*temp;
		return rV;
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

/////////////////////////////////////////////////////////////////////////////
class normalLLTESTA3 : public normalLL {

public:
	normalLLTESTA3() {};
	normalLLTESTA3(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
		bool SS) : normalLL(tY, tX, SS) {
		// if it is a sufficient statistics model
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

	// mean and variance for every dose group
	int    nParms() {
		return meanX.cols()+2;
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		Eigen::MatrixXd temp(theta.size() - 2, 1);
		for (int i = 0; i < temp.size(); i++) {
			temp(i, 0) = theta(i, 0);
		}
		 
		Eigen::MatrixXd rV = meanX * temp;
		return rV;
	};

	Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return variance(theta, X);
	};

	Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		Eigen::MatrixXd tmean = mean(theta,d);
		Eigen::MatrixXd rV = exp(theta(theta.rows() - 1, 0))*pow(abs(tmean.array()), theta(theta.rows() - 2, 0));//  theta(theta.rows() - 2, 0));
		return rV;
	};

public:
	std::vector<double> udoses; // unique dose groups 
	Eigen::MatrixXd     meanX;

};

/////////////////////////////////////////////////////////////////////////////
class normalLLTESTR : public normalLL {

public:
	normalLLTESTR() {};
	normalLLTESTR(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
		bool SS) : normalLL(tY, tX, SS) {
	
	};

	// mean and variance for every dose group
	int    nParms() {
		return 2;
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta) {
		return mean(theta, X);
	};

	virtual Eigen::MatrixXd mean(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		Eigen::MatrixXd rV = d; 
		rV = rV.array()*0 +  theta(0, 0);
		return rV;
	};

	Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
		return variance(theta, X);
	};

	Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		Eigen::MatrixXd tmean = mean(theta, d);
		Eigen::MatrixXd rV = tmean.array() * 0 + exp(theta(1, 0)); 
		return rV;
	};

public:
	std::vector<double> udoses; // unique dose groups 
	Eigen::MatrixXd     meanX;

};

#endif
