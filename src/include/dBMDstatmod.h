//File: dBMDstatmod.h
//Purpose: Creates a basic statistical model class
//         that is a BMD statistical model. This class provides
//         all of the information necessary to compute a benchmark dose
//         for dichotomous data. 
//Creator: Matt Wheeler
//Date   : 12/21/2017
//Changes: 
//       Date: 11/06/2018 - Added code in the profile function to stop 
//							the upper bound calculation when it is 5 times greater
//					        than the maximum administered dose. 
//

#pragma once

#ifndef dBMDstatmodH
#define dBMDstatmodH

#include "statmod.h"
#include <cmath>
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

#include <list>
//#include <cmath>
#include <nlopt.hpp>
#include <limits>
#include <iostream>

using namespace std; 
//using namespace Rcpp; // Only needed for debugging in R




template <class LL, class PR>
class dBMDModel : public statModel<LL,PR> {
public:
	dBMDModel(LL t_L, PR t_PR,
			 std::vector<bool> b_fixed,
			 std::vector<double> d_fixed):statModel<LL,PR>(t_L,t_PR,b_fixed,d_fixed){
	};

	std::vector<double > fixConstraintLB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double lb) {
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.fixConstraintLB(theta, BMD, BMR, isExtra, lb);
	}

	std::vector<double > fixConstraintUB(Eigen::MatrixXd theta, double BMD, double BMR, bool isExtra, double ub) {
		
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.fixConstraintUB(theta, BMD, BMR, isExtra, ub);
	}

	double BMR_CONSTRAINT(Eigen::MatrixXd theta, double *grad, double BMR, bool isExtra) {
		
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.BMR_CONSTRAINT(theta, grad, BMR,isExtra);
	}

	// These functions return the added risk and extra risk
	// NC is a no covariates option 
	double added_riskBMDNC(double BMR) {
		// set all the fixed parameters to their value
		// to do check if the model has an estimate
		return statModel<LL, PR>::log_likelihood.compute_BMD_ADDED_NC(statModel<LL, PR>::theta,BMR);
	};

	double extra_riskBMDNC(double BMR) {
		// to do check if the model has an estimate
		// set all the fixed parameters to their value
		return statModel<LL, PR>::log_likelihood.compute_BMD_EXTRA_NC(this->getEST(),BMR);
	};

	Eigen::MatrixXd start_equalityOptim(Eigen::MatrixXd theta, double BMD, double BMR,bool isExtra) {
		
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		
		if (isExtra) {
			return statModel<LL, PR>::log_likelihood.beta_BMD_ExtraNC(theta, BMD, BMR);
		}
		else {
			return statModel<LL, PR>::log_likelihood.beta_BMD_AddedNC(theta, BMD, BMR);
		}
	}

	Eigen::MatrixXd bmdAddedNC_equals(Eigen::MatrixXd theta,double BMD,double BMR) {
		return statModel<LL, PR>::log_likelihood.beta_BMD_AddedNC(theta, BMD, BMR);
	}

	Eigen::MatrixXd bmdExtraNC_equals(Eigen::MatrixXd theta, double BMD, double BMR) {
		return statModel<LL, PR>::log_likelihood.beta_BMD_ExtraNC(theta, BMD, BMR);
	}

	bool fixedNC_BMDformula() {
		return statModel<LL, PR>::log_likelihood.fixedNC_BMDformula();
	}

	int removedParameter() {
		return statModel<LL, PR>::log_likelihood.parameterRemoved();
	}

	double equality_extra(Eigen::MatrixXd theta, double *grad, double BMD, double BMR) {
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.compute_BMD_EXTRA_NC_EQUALITY(theta, grad,BMD, BMR);
	}
	double equality_added(Eigen::MatrixXd theta, double *grad, double BMD, double BMR) {
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.compute_BMD_ADDED_NC_EQUALITY(theta,grad, BMD, BMR);
	}

	virtual double inequality_extra(Eigen::MatrixXd theta, double BMD, double BMR, double inequality,bool geq, double *grad) {
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.compute_BMD_EXTRA_NC_INEQUALITY(theta, BMD, BMR,inequality,geq, grad);
	}

	virtual double inequality_added(Eigen::MatrixXd theta, double BMD, double BMR, double inequality, bool geq,double * grad) {
		// set all the fixed parameters to their value
		for (int i = 0; i < this->isFixed.size(); i++) {
			if (this->isFixed[i]) {
				theta(i, 0) = this->fixedV[i];
			}
		}
		return statModel<LL, PR>::log_likelihood.compute_BMD_ADDED_NC_INEQUALITY(theta, BMD, BMR,inequality,geq,grad);
	}


};

template <class A, class B>
struct optimInfo {
	dBMDModel<A, B> * sm;
	double cBMD;
	double BMR;
	bool   isExtra;
};

template <class A, class B>
struct inequalityInfo : public optimInfo<A, B> {

	double inequality; 
	bool geq; 
};

/////////////////////////////////////////////////////////////////
//Function: DICHOTOMOUS_BMD_neg_pen_likelihood(unsigned n,
//								const double *b,
//								double *grad,
//								void *data)
//Purpose: Used to estimate the negative penalized likelihood from the statMod class
//         in the proper NLOPT format
//Input: 
//         unsigned n: The length of the vector of parameters to the optimizer
//         const double *b   : A pointer to the place in memory the vector resides
//         double *gradient  : The gradient of the function at b, NULL if unused. 
//                             It is currently unused. 
//         void    *data     : Extra data needed. In this case, it is a statModel<LL,PR> object,
//							   which is used to compute the negative penalized likelihood
//////////////////////////////////////////////////////////////////
template <class LL, class PR>
double DICHOTOMOUS_BMD_neg_pen_likelihood(unsigned n,
						   const double *b,
						   double *grad,
						   void *data) {
	optimInfo<LL,PR> *model = (optimInfo<LL,PR> *) data;
	Eigen::MatrixXd theta(n, 1);
	for (int i = 0; i < n; i++) theta(i, 0) = b[i];
		
	if (model->isExtra) {
		theta = model->sm->bmdExtraNC_equals(theta, model->cBMD, model->BMR);
	}
	else {
		theta = model->sm->bmdAddedNC_equals(theta, model->cBMD, model->BMR);
	} 
	// if an optimizer is called that needs a gradient we calculate
	// it using the finite difference gradient
	if (grad) {
		Eigen::MatrixXd mgrad = model->sm->gradient(theta);

		int rp = model->sm->removedParameter();
		// Copy over the bounds minus the "extra parameter"
		int j = 0;
		for (int i = 0; i < model->sm->nParms(); i++) {
			if (i != rp) {
				grad[j] = mgrad(i, 0);
				j++;
			}
		}
	}

	double rV = model->sm->negPenLike(theta);
//cout << rV << endl; 
	return  rV;
}


/////////////////////////////////////////////////////////////////
//Function: equality_constraint(unsigned n,
//								const double *b,
//								double *grad,
//								void *data)
//Purpose: Used to estimate the negative penalized likelihood from the statMod class
//         in the proper NLOPT format
//Input: 
//         unsigned n: The length of the vector of parameters to the optimizer
//         const double *b   : A pointer to the place in memory the vector resides
//         double *gradient  : The gradient of the function at b, NULL if unused. 
//                             It is currently unused. 
//         void    *data     : Extra data needed. In this case, it is a statModel<LL,PR> object,
//							   which is used to compute the negative penalized likelihood
//////////////////////////////////////////////////////////////////
template <class LL, class PR>
double equality_constraint(unsigned n,
	const double *b,
	double *grad,
	void *data) {
	optimInfo<LL, PR> *model = (optimInfo<LL, PR> *) data;
	Eigen::MatrixXd theta(n, 1);
	for (int i = 0; i < n; i++) theta(i, 0) = b[i]; // simple addition

	if (model->isExtra) {
		return  model->sm->equality_extra(theta,grad,model->cBMD,model->BMR);
	}
	else {
		return  model->sm->equality_added(theta,grad,model->cBMD, model->BMR);
	}

	
}

template <class LL, class PR>
double inequality_constraint(unsigned n,
	const double *b,
	double *grad,
	void *data)
{
	inequalityInfo<LL,PR> *model = (inequalityInfo<LL, PR> *) data;
	Eigen::MatrixXd theta(n, 1);
	for (int i = 0; i < n; i++) theta(i, 0) = b[i]; 

	if (model->isExtra) {
		return  model->sm->inequality_extra(theta, model->cBMD, model->BMR,model->inequality,model->geq,grad);
	}else {
		return  model->sm->inequality_added(theta, model->cBMD, model->BMR, model->inequality, model->geq,grad);
	}


}


////////////////////////////////////////////////////////////
// Function: inequality_constraint_general(unsigned n,
//											const double *b,
//											double *grad,
//										    void *data)
//Purpose: Many BMR calculations require functions of parameters to be greater than zero
//         this sets up this requirement (e.g. log(f(k)) where f(k) is a function of a vector
//         of k parameters). 
//Input  : unsigned n: Number of parameters
//         const double *b : pointer to the parameter values
//		   double *grad    : gradient 
//         void * data     : additional data needed for the inequality
//Output : Calls the function M->BMR_CONSTRAINT
template <class LL, class PR>
double inequality_constraint_general(unsigned n,
				     const double *b,
				     double *grad,
				     void *data)
{
	inequalityInfo<LL, PR> *model = (inequalityInfo<LL, PR> *) data;
	Eigen::MatrixXd theta(n, 1);
	for (int i = 0; i < n; i++) theta(i, 0) = b[i]; // simple addition

	if (model->isExtra) {
		theta = model->sm->bmdExtraNC_equals(theta, model->cBMD, model->BMR);
	}
	else {
		theta = model->sm->bmdAddedNC_equals(theta, model->cBMD, model->BMR);
	}

	return   model->sm->BMR_CONSTRAINT(theta, grad, model->BMR, model->isExtra);

}




template <class LL, class PR>
optimizationResult findMAX_W_EQUALITY( dBMDModel<LL, PR>  *M,
				       Eigen::MatrixXd start,
				       const bool isExtra, // true if it is false if it is added
				       double BMD,
				       double BMR) {
	
	optimizationResult oR;
	std::vector<double> vec_start(start.rows()); 
	for (int i = 0; i < vec_start.size(); i++) {
		vec_start[i] = start(i, 0); 
	}

	Eigen::MatrixXd temp_data = M->parmLB();
	std::vector<double> lb(M->nParms());
	for (int i = 0; i < M->nParms(); i++) lb[i] = temp_data(i, 0);
	temp_data = M->parmUB();
	std::vector<double> ub(M->nParms());
	for (int i = 0; i < M->nParms(); i++) ub[i] = temp_data(i, 0);

	// you have to start on the equality constraint
	std::vector<double> x(M->nParms());
	

	start = M->start_equalityOptim(start, BMD, BMR, isExtra);
	for (int i = 0; i < M->nParms(); i++) x[i] = start(i, 0); 

	/// check the box-bounds satisfy the equality constraint


	//This is all very specific to the multistage
	//which makes this code not generalizable to models requiring equality constraints
	// however, the multistage model is the only code that runs 
	// this code... fix this in the future for generalizability
//	for (double g:x){cout << g << endl; }
	bool lbsatisfied = true; bool ubsatisfied = true; 
	for (int i = 0; i < x.size(); i++) {
		if (x[i] < lb[i]) { lbsatisfied = false; }
		if (x[i] > ub[i]) { ubsatisfied = false; }
	}
	double g = 1.0 / (1.0 + exp(-x[0])); 
	if (!lbsatisfied || !ubsatisfied) {
		x = vec_start;   // don't worry about the gamma parameter 
      
		
		double t_sum = 0.0; 
		if (isExtra) {
			 // worry about the 
			for (int j = 1; j < x.size()-1; j++) {
				
				t_sum += x[j] * pow(BMD, j);
				
			}
			x[x.size()-1] = ((-log(1.0 - BMR) - t_sum) / pow(BMD, x.size()-1));
		}
		else {
			for (int j = 1; j < x.size()-1; j++) {
				t_sum += x[j] * pow(BMD, j);	
			}
			x[x.size()-1] = ((-log(1.0 - BMR / (1 - g)) - t_sum) / pow(BMD, x.size()-1));
		}

	}
	
	lbsatisfied = true; ubsatisfied = true; 
	for (int i = 0; i < x.size(); i++) {
		if (x[i] < lb[i]) { lbsatisfied = false; }
		if (x[i] > ub[i]) { ubsatisfied = false; }
	}
	
	if (!lbsatisfied || !ubsatisfied){
		for (int i =1; i < x.size(); i++){
			x[i] = 0.0; 
		}
		if (isExtra){
			x[1] = (-log(1.0 - BMR)) / BMD; 
		} else{
			x[1] = (-log(1.0 - BMR / (1 - g))) / BMD; 
		}
	}
	

	double minf;
			
	nlopt::opt opt(nlopt::LN_AUGLAG, M->nParms());
	
	nlopt::opt local_opt(nlopt::LD_LBFGS,M->nParms()); // set up the local optimizer
	nlopt::opt local_opt2(nlopt::LN_SBPLX,M->nParms()); 


	local_opt.set_xtol_abs(1e-3);
	local_opt.set_initial_step(1e-4); 
	local_opt.set_maxeval(10000); 
	
	local_opt2.set_xtol_abs(1e-3);
	local_opt2.set_initial_step(1e-4); 
	local_opt2.set_maxeval(10000); 
	
	local_opt.set_lower_bounds(lb); 
	local_opt2.set_lower_bounds(lb); 
	local_opt.set_upper_bounds(ub); 
	local_opt2.set_upper_bounds(ub); 

	bool good_opt = false; 
	int opt_iter = 0;
	///////////////////////////////////
	///////////////////////////////////////////////////////////////
	// create the information for the local equality constraint
	optimInfo<LL, PR> info;
	info.BMR = BMR;
	info.cBMD = BMD;
	info.sm = M;
	info.isExtra = isExtra;
	
	////////////////////////////////////////////////
	//////////////////////////////////////////////
	// set the start distance to a size
	// nlopt's default options suck
	std::vector<double> init(x.size());
	for (int i = 0; i < x.size(); i++) init[i] = 1e-4;
	opt.set_initial_step(init);
	local_opt.set_initial_step(init);

	///////////////////////////////////////////////
	nlopt::result result = nlopt::FAILURE; 
	opt.add_equality_constraint(equality_constraint <LL, PR>, &info, 1e-4);
	opt.set_min_objective(neg_pen_likelihood<LL, PR>, M);
    //ofstream file;
    //file.open("bmds.log", fstream::app);
    //file << __FUNCTION__ << " at line: " << __LINE__ << " Before optimize"<< endl;
	while( opt_iter < 2 &&  !good_opt){
		opt_iter++; 
		if (opt_iter == 0)
			opt.set_local_optimizer((const nlopt::opt) local_opt);
		else
			opt.set_local_optimizer((const nlopt::opt) local_opt2);

		opt.set_lower_bounds(lb); 
		opt.set_upper_bounds(ub);
		opt.set_xtol_abs(1e-4);
		opt.set_maxeval(20000);
        // Ensure that starting values are within bounds
        for (int i = 0; i < x.size(); i++) {
          double temp = x[i];
          if (temp < lb[i]) temp = lb[i];
          else if (temp > ub[i]) temp = ub[i];
          x[i] = temp;
        } // end for

		try {
			result = opt.optimize(x, minf);
			good_opt = true; //optimization succeded
		}catch (nlopt::roundoff_limited) {
			good_opt = false;
//cerr << "Round Off Limited" << endl;  
		}catch (nlopt::forced_stop) {
			good_opt = false;
//cerr << "Forced Stop" << endl; 
		}
		catch (const std::invalid_argument &exc) {
          //file << "\tline " << __LINE__ << ": invalid arg, opt_iter= " << opt_iter << endl;
          //flush(file);
//cerr << "Invalid Argument" << endl; 
			good_opt = false;
		}catch (const std::exception &exc) {
			good_opt = false;
//cerr << "Std Exeption" << endl; 
		}catch (...) { 
			good_opt = false;
//cerr << "default exception" << endl; 
		}

		if (result > 5) { // Either 5 = 
			good_opt = false;
		}
		
        //file.close();

	}
	
	Eigen::Map<Eigen::MatrixXd> d(x.data(), M->nParms(), 1); // return values
	oR.result = result;
	oR.functionV = minf;
	oR.max_parms = d;
	return oR;
}


////////////////////////////////////////////////////////////////////
//Function: DICHOTOMOUS_findMAX(statModel<LL, PR>  *M, 
//						Eigen::MatrixXd start,
//						double BMD,
//						double BMR)
//Find the maximum a posteriori estimate given a statistical 
//model
//Input : statMod<LL,PR> *M - A given dichotomous statistical model
//                          - start: A starting value
//                          - BMD  : The BMD target
//                          - BMR  : The BMR
//Output: statMod<LL,PR> *M - The model with it's MAP parameter set.

////////////////////////////////////////////////////////////////////
//Function: DICHOTOMOUS_findMAX(statModel<LL, PR>  *M, 
//						Eigen::MatrixXd start,
//						double BMD,
//						double BMR)
//Find the maximum a posteriori estimate given a statistical 
//model
//Input : statMod<LL,PR> *M - A given dichotomous statistical model
//                          - start: A starting value
//                          - BMD  : The BMD target
//                          - BMR  : The BMR
//Output: statMod<LL,PR> *M - The model with it's MAP parameter set.
template <class LL, class PR>
optimizationResult findMAX(    dBMDModel<LL, PR>  *M,
								Eigen::MatrixXd start, 
								const bool isExtra, // true if it is false if it is added
								double BMD,
								double BMR,
								const int iter, 
								nlopt::algorithm algorithm) {
	optimizationResult oR;
	if (!M->fixedNC_BMDformula()) {
		// there is no way to express one of the parameters as a function of the others
		// thus we have to do this using equality constraints
		return findMAX_W_EQUALITY<LL,PR>(M, start, isExtra, BMD, BMR);
	}
		
		optimInfo<LL, PR> info;
		info.BMR = BMR;
		info.cBMD = BMD;
		info.sm = M;
		info.isExtra = isExtra;
		inequalityInfo<LL, PR> ineqinfoLB, ineqinfoUB;
		ineqinfoLB.BMR = BMR;
		ineqinfoLB.cBMD = BMD;
		ineqinfoLB.sm = M;
		ineqinfoLB.isExtra = isExtra; 

		ineqinfoUB.BMR = BMR;
		ineqinfoUB.cBMD = BMD;
		ineqinfoUB.sm = M;
		ineqinfoUB.isExtra = isExtra;

		Eigen::MatrixXd temp_data = M->parmLB();
		int rp = M->removedParameter();
		std::vector<double> lb(M->nParms() - 1);

		nlopt::opt opt(algorithm, M->nParms() - 1);

		// Copy over the bounds minus the "extra parameter"
		int j = 0;
		for (int i = 0; i < M->nParms(); i++) {
			if (i != rp) {
				lb[j] = temp_data(i, 0);
				j++;
			}
			else {
				ineqinfoLB.geq = true;
				ineqinfoLB.inequality = temp_data(i, 0);
				opt.add_inequality_constraint(inequality_constraint <LL, PR>, &ineqinfoLB );
				
			}
		}

		temp_data = M->parmUB();
		std::vector<double> ub(M->nParms() - 1);

		// Copy over the bounds minus the "extra parameter"
		j = 0;
		
		for (int i = 0; i < M->nParms(); i++) {
			if (i != rp) {
				ub[j] = temp_data(i, 0);
				j++;
			}
			else {
				ineqinfoUB.geq = false;
				ineqinfoUB.inequality = temp_data(i, 0);
				opt.add_inequality_constraint(inequality_constraint <LL, PR>, &ineqinfoUB);
			}

		}
		// check inequality constraints	

		std::vector<double> x(M->nParms() - 1);
		if (start.rows() == M->nParms()) { // we need to remove the extra parameter
			int j = 0;
			for (int i = 0; i < M->nParms(); i++) {
				if (i != rp) {
					x[j] = start(i, 0);
					j++;
				}
				
			}
		}
		else { // our starting values are assumed to have the extra parameter removed
			for (int i = 0; i < start.rows(); i++) {
				x[i] = start(i, 0);
			}
		}


		Eigen::MatrixXd temp(x.size(), 1); 
		for (int i = 0; i < temp.rows(); i++) temp(i, 0) = x[i]; 


		////////////////////////////////////////////////////////////////////////////
		// Check to see if our starting value has the inequality constraints violated
		double t1, t2; 
		if (isExtra) {
			t1 = M->inequality_extra(temp, BMD, BMR, ineqinfoLB.inequality, true, NULL); 
			t2 = M->inequality_extra(temp, BMD, BMR, ineqinfoUB.inequality, false, NULL);
		}
		else {
			t1 = M->inequality_added(temp, BMD, BMR, ineqinfoLB.inequality, true, NULL);
			t2 = M->inequality_added(temp, BMD, BMR, ineqinfoUB.inequality, false, NULL);
		}
		//for (double g:x){ cout << g << endl; }
			
		if (t1 > 0.0 || t2 > 0.0) {
			
			//constraints are violated need to do something
			if (t1 > 0.0) {
				x = M->fixConstraintLB(temp, BMD, BMR, isExtra, ineqinfoLB.inequality);
			}
			else {
				x = M->fixConstraintUB(temp, BMD, BMR, isExtra, ineqinfoUB.inequality);
			}
		}
		//////////////////////////////////////////////////////////////////
		//for (double g:x){ cout << g << endl; }
		//cout << BMD << endl; 

		opt.add_inequality_constraint(inequality_constraint_general<LL, PR>, &ineqinfoUB);

		double minf;
		
		opt.set_lower_bounds(lb);
		opt.set_upper_bounds(ub);
		opt.set_ftol_rel(1e-3);

		opt.set_maxeval(iter);

		opt.set_min_objective(DICHOTOMOUS_BMD_neg_pen_likelihood<LL, PR>, &info);

		//////////////////////////////////////////////
		// set the start distance to a size
		// nlopt's default options suck
		std::vector<double> init(x.size());
		for(int i = 0; i < x.size(); i++) init[i] = 1e-4; 
		opt.set_initial_step(init);
		///////////////////////////////////////////////
	
		nlopt::result result = nlopt::FAILURE;
		result =  opt.optimize(x, minf);
		
		Eigen::Map<Eigen::MatrixXd> d(x.data(), M->nParms() - 1, 1);
		oR.result = result;
		oR.functionV = minf;
		oR.max_parms = d;
		return oR;
}

///////////////////////////////////////////////////////////////////////////////
//Function fit_profileLogic(dBMDModel<LL, PR>  *M,
//						 bool isExtra,		// true if it is false if it is added
//						 double BMR,
//		                 double BMDchange,
//						 double totalChange)
//Purpose: This function fits a the profile likelihood using standard logic
//Input  : dBMDModel<LL, PR>  *M - Dichotomous BMD model   
//		   bool isExtra          - true if it is extra risk	
//         double BMD			 - current BMD at the MAP
//		   double BMR            - BMR [0,1]
//	       double BMDchange      - %Change in the BMD to increment each time
//		   double totalChange    - totalChange in penalized likelihood before one stops
// Output: Returns a list of vectors.  The first item is the fit information.
//                                     The second item is the parameter optimization values.
//									   The third item is the list of parameter optimization values including the 
//                                     fixed parameter value. 
//                                     On an error there is only one item in the list indicating the return 
//                                     code. 
////////////////////////////////////////////////////////////////////////////////////
template <class LL, class PR>
std::list<Eigen::MatrixXd> fit_profileLogic(dBMDModel<LL, PR>  *M,
											 Eigen::MatrixXd start,
											 const bool isExtra, // true if it is false if it is added
											 double BMD,
										     double BMR,
											 const int iter,
											 nlopt::algorithm algorithm) {
	optimizationResult oR;
	std::list<Eigen::MatrixXd> rV; 
	Eigen::MatrixXd result(3, 1); 
	Eigen::MatrixXd mParms, cParms; //mParms - Maximized parms in the restricted model
							        //cParms - Complete  model parameters including the fixed value
	bool bad_result = false;
	
	// try to fit the model
	try {
	
//	cout << start << endl << "BMD:" << BMD << endl; 
		oR = findMAX<LL, PR>(M,
							start,
							isExtra,
							BMD, BMR, iter,
							algorithm);
	
	}
	catch (nlopt::roundoff_limited) {
		//TODO check to see if their is a roundoff limitation we 
		// mark the failure and set the maximum function value to infinity. 
		result(0, 0) = DBL_MAX; result(1, 0) = BMD; result(2, 0) = NLOPT_ROUNDOFF_LIMITED;
		bad_result = true;
		rV.push_front(result);
	//	cout << "AB" << endl;
		return rV;
	}
	catch (nlopt::forced_stop) {
		result(0, 0) = DBL_MAX; result(1, 0) = BMD; result(2, 0) = NLOPT_FORCED_STOP;
		bad_result = true;
		rV.push_front(result);
	//	cout << "BB" << endl;
		return rV;
	}
	catch (const std::invalid_argument &exc) {
		result(0, 0) = DBL_MAX; result(1, 0) = BMD; result(2, 0) = NLOPT_FAILURE;
		bad_result = true;
		rV.push_front(result);
	//	cout <<"CB" << endl;
		return rV;
	}
	catch (const std::exception &exc) {
        	result(0, 0) = DBL_MAX; result(1, 0) = BMD; result(2, 0) = NLOPT_FAILURE;
		bad_result = true;
		rV.push_front(result);
	//	cout <<"DB" << endl;
		return rV;
        }catch (...){
		result(0, 0) = DBL_MAX; result(1, 0) = BMD; result(2, 0) = NLOPT_FAILURE;
		bad_result = true;
		rV.push_front(result);
	//	cout <<"CC" << endl;
		return rV;
	}

	// if there was no problem in the optimization
	if (!bad_result) {
		result(0, 0) = oR.functionV; result(1, 0) = BMD; result(2, 0) = oR.result;
		mParms = oR.max_parms; 
		if (isExtra) {
			cParms = M->bmdExtraNC_equals(mParms, BMD, BMR);
		}
		else {
			cParms = M->bmdAddedNC_equals(mParms, BMD, BMR);
		}
		rV.push_front(cParms);
		rV.push_front(mParms);
	}
	
	rV.push_front(result);
	return rV; 
}

///////////////////////////////////////////////////////////////////////////////
//Function profile_BMDNC(dBMDModel<LL, PR>  *M,
//						 bool isExtra,		// true if it is false if it is added
//						 double BMR,
//		                 double BMDchange,
//						 double totalChange,
//						 bool robust)
//Purpose: This function iteratively changes the BMD by a BMDchange%  
//         until a total change in the penalized likelihood is found. 
//Input  : dBMDModel<LL, PR>  *M - Dichotomous BMD model   
//		   bool isExtra          - true if it is extra risk	
//         double BMD			 - current BMD at the MAP
//		   double BMR            - BMR [0,1]
//	       double BMDchange      - %Change in the BMD to increment each time
//		   double totalChange    - totalChange in penalized likelihood before one stops
//         bool   robust         - true if we do a robust search of the optimization space, false otherwise
////////////////////////////////////////////////////////////////////////////////
template <class LL, class PR>
Eigen::MatrixXd profile_BMDNC(  dBMDModel<LL, PR>  *M,
								const bool isExtra,		// true if it is false if it is added
								const double BMD,
								const double BMR,
								const double BMDchange,
								const double totalChange, 
								bool robust)
								
{
	//double mapBMD = BMD; // current BMD evaluated at estimate M->getEST 
	Eigen::MatrixXd parms = M->getEST();
	Eigen::MatrixXd X = M->returnX(); 
	

	double max_dose = X.maxCoeff(); 

	double MAP_LIKE = M->negPenLike(parms);
	Eigen::MatrixXd ret1(3,1), ret2, fParms; 
	std::list<Eigen::MatrixXd> CL,temp1,temp2; 

	double CLIKE = MAP_LIKE; 	
	int max_iters = 0; 
	double CBMD = BMD * (1.0-BMDchange);
	double PLIKE = CLIKE;
	
	bool error = false;
	
	ret1(0, 0) = PLIKE; ret1(1, 0) = BMD; ret1(2, 0) = 666; 
	CL.push_front(ret1); 
	//
	// Start off going down until the criteria is satified
	while (fabs(MAP_LIKE - CLIKE) < totalChange && max_iters < 500 && CBMD > 1e-8)
	{
		
		// fit the first profile likelihood option			
		temp1 = fit_profileLogic<LL, PR>(M,
			parms,
			isExtra,
			CBMD, BMR, 10000, nlopt::LN_COBYLA);
	       
		// if the profile fit is called using the 
		// robust option then we call the second optimizer 
		// second opinion
		if (true) {
			temp2 = fit_profileLogic<LL, PR>(M,
				parms,
				isExtra,
				CBMD, BMR, 250, nlopt::LD_MMA);
			ret1 = temp1.front();
			ret2 = temp2.front();

			if (ret2(0, 0) < ret1(0, 0)) // second optimization was "the best"
			{
				 
				if (temp2.size() == 1) { // there was an error in the optimization 
					error = true;
					for (Eigen::MatrixXd n : temp1) {
						ret1 = n;
					}
				}
				else {
					int k = 0;
					for (Eigen::MatrixXd n:temp2) {
						switch (k) {
							case 1: parms = n; 	 break;
							case 2: fParms = n;  break;
							default: ret1 = n;   break;
						}
						k++;
					}
				}
			}
			else { // first optimization was the best
			       
				if (temp1.size() == 1) { // there was an error in the optimization 
					error = true;
					
					for (Eigen::MatrixXd n : temp1) {
						ret1 = n;
					}
				}
				else {
										
					int k = 0;
					for (Eigen::MatrixXd n:temp1) {
						switch (k) {
							case 1: parms = n; 	 break;
							case 2: fParms = n;  break;
							default: ret1 = n;   break;
						}
						k++;
					}
				}
			}

		}
		else {
			if (temp1.size() == 1) { // there was an error in the optimization 
				error = true;
				for (Eigen::MatrixXd n : temp1) {
					ret1 = n;
				}
			}
			else {
				int k = 0;
				for (Eigen::MatrixXd n:temp1) {
					switch (k) {
					case 1: parms = n; 	 break;
					case 2: fParms = n;  break;
					default: ret1 = n;   break;
					}
					k++;
				}
			}

		}
		
	//	
		CLIKE = ret1(0,0); 
		PLIKE = CLIKE;
		CBMD *= 1-BMDchange;
		CL.push_front(ret1);
		max_iters++; 
		if (error) { // there was an error in the optimization we break
			break; 
		}
	}
	
	

	CLIKE = MAP_LIKE;  	CBMD = BMD * (1.0 + BMDchange);
	int max_iters2 = 0;   	parms = M->getEST();
        error = false; 
	
	// Now increase the BMD sequentially going up
 // stop when we are 2.5 times greater than the maximum dose
	while ( max_iters2 < 200
	  	 && CBMD < 2.5*max_dose
	  	 && fabs(MAP_LIKE - CLIKE) < totalChange )
	{   
		// fit the first profile likelihood option			
		temp1 = fit_profileLogic<LL, PR>(M,
						parms,
						isExtra,
						CBMD, BMR, 10000, nlopt::LN_COBYLA);
			//			CBMD, BMR, 10000,         nlopt::LD_MMA);
		// if the profile fit is called using the 
		// robust option then we call the second optimizer 
		// second opinion
		if (robust) {
			temp2 = fit_profileLogic<LL, PR>(M,
				parms,
				isExtra,
				CBMD, BMR, 150, nlopt::LD_MMA);

			ret1 = temp1.front();
			ret2 = temp2.front();

			if (ret2(0, 0) < ret1(0, 0)) // second optimization was "the best"
			{
				if (temp2.size() == 1) { // there was an error in the optimization 
					error = true;
					for (Eigen::MatrixXd n : temp1) {
						ret1 = n;
					}
				}
				else {
					int k = 0;
					for (Eigen::MatrixXd n:temp2) {
						switch (k) {
						case 1: parms = n; 	 break;
						case 2: fParms = n;  break;
						default: ret1 = n;   break;
						}
						k++;
					}
				}
			}
			else { // first optimization was the best
				
				if (temp1.size() == 1) { // there was an error in the optimization 
					error = true;
					for (Eigen::MatrixXd n : temp1) {
						ret1 = n;
					}
				}
				else {
					int k = 0;
					for (Eigen::MatrixXd n:temp1) {
						switch (k) {
						case 1: parms = n; 	 break;
						case 2: fParms = n;  break;
						default: ret1 = n;   break;
						}
						k++;
					}
				}
			}

		}
		else {
			if (temp1.size() == 1) { // there was an error in the optimization 
				error = true;
				for (Eigen::MatrixXd n : temp1) {
					 ret1 = n;   
				}
			}
			else {
				int k = 0;
				for (Eigen::MatrixXd n:temp1) {
					switch (k) {
					case 1: parms = n; 	 break;
					case 2: fParms = n;  break;
					default: ret1 = n;   break;
					}
					k++;
				}
			}

		}
       
               
		CLIKE = ret1(0, 0);
		PLIKE = CLIKE;
		CBMD *= 1 + BMDchange;
	//	cout << "ARG:" << CBMD << endl; 
		if (error) { // there was an error in the optimization we break
			break;
		}else{
			CL.push_back(ret1);
		}
		
	
		max_iters2++;
	}	
	
	// Take the lists and put them in one big
	// matrix to return. 
	Eigen::MatrixXd returnMat(CL.size(), 3);
	int ii = 0; 
	for (Eigen::MatrixXd it:CL) {
		returnMat.row(ii) = it.transpose(); 
		ii++; 
	}
	//Rcpp::Rcout << returnMat << endl; 
	returnMat.col(0).array() =  round(round(1e4*returnMat.col(0).array()) - round(1e4*MAP_LIKE))/1e4; // round it to be within the optimizers tol
	
	return returnMat; 
}



#endif

