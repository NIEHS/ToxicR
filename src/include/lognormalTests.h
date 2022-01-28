#pragma once

#ifndef lognormalTestsH
#define lognormalTestsH

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

#include <gsl/gsl_randist.h>
#include <log_likelihoods.h>
#include <lognormal_likelihoods.h>

/*implements the basic model model Y = X*/
class lognormalLLTESTA1 : public lognormalLL {

			public:
				lognormalLLTESTA1() {};
				lognormalLLTESTA1(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
												   bool SS) : lognormalLL(tY, tX, SS){
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
					return log(rV.array());
				}; 

				Eigen::MatrixXd variance(Eigen::MatrixXd theta) {
					return variance(theta, X);
				}; 

				 Eigen::MatrixXd variance(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
					double var = exp(theta(theta.rows() - 1, 0));
					Eigen::MatrixXd rV = Eigen::MatrixXd::Ones(d.rows(), 1)*var;
					return rV;
				}; 


			public:
					std::vector<double> udoses; // unique dose groups 
					Eigen::MatrixXd     meanX; 

};

////////////////////////////////////////////////
class lognormalLLTESTA2 : public lognormalLL {

public:
	lognormalLLTESTA2() {};
	lognormalLLTESTA2(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
		bool SS) : lognormalLL(tY, tX, SS) {
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
		return log(rV.array());
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

public:
	std::vector<double> udoses; // unique dose groups 
	Eigen::MatrixXd     meanX;

};


/////////////////////////////////////////////////////////////////////////////
class lognormalLLTESTR : public lognormalLL {

public:
  lognormalLLTESTR() {};
  lognormalLLTESTR(Eigen::MatrixXd tY, Eigen::MatrixXd tX,
    bool SS) : lognormalLL(tY, tX, SS) {

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
    rV = rV.array() * 0 + theta(0, 0);
    return log(rV.array());
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
