//File: statisticalmodel.h
//Purpose: Create a basic statistical model class
//		   that is used for any general analysis, and can be 
//         rather flexible. 
//Creator: Matt Wheeler
//Date   : 12/18/2017
//Changes: 
//
//
//
#pragma once
#ifndef statisticalmodelH
#define statisticalmodelH

#include <nlopt.hpp>
#ifdef R_COMPILATION
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
   
#endif

#include <gsl/gsl_randist.h>

template <class LL, class PR>
class statModel {
public:
	//basic constructor
	statModel(LL t_L, PR t_PR) :log_likelihood(t_L), prior_model(t_PR) {

	};

private:
	// A stat model has a Log Likelihood and 
	// A prior model over the parameters
	// Note: The prior over the parameters does not have to be an actual 
	// prior. It can merely place bounds (e.g., box, equality) on functions
	// of the parameters. 
	LL log_likelihood;
	PR prior_model;

};
#endif

