#include "stdafx.h" // Precompiled header - does nothing if building R version
#include "lognormal_HILL_NC.h"

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include "gsl/gsl_cdf.h"
#include <gsl/gsl_randist.h>
    
#include <iostream>
#include <cmath>
#include <math.h>      

using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////
// function: parameter_to_remove()
// purpose: Tell the optimizer which profile likelihood method is
//			best for the given bmdtype.  For all models, the HYBRID
//          is always represented as an equality constraint. 
// input: 
//	contbmd TYPE
// output: 
//  PROFILE_INEQUALITY - One of the parameters can be made equal to the others as a function
//						 of the fixed BMD. The optimizer thus optimizes a smaller problem
//  PROFILE_EQUALITY   - The BMD is  a function of multiple parameters and can not be disentangled
//                       An equality constraint is used here. 
///////////////////////////////////////////////////////////////////////////////////////////////////
int lognormalHILL_BMD_NC::parameter_to_remove(contbmd TYPE) {

	switch (TYPE) {
	case CONTINUOUS_BMD_ABSOLUTE:
		return 1; // slope
	case CONTINUOUS_BMD_STD_DEV:
		return nParms() - 1;
		break;
	case CONTINUOUS_BMD_REL_DEV:
	case CONTINUOUS_BMD_POINT:
	case CONTINUOUS_BMD_EXTRA:
		return 0;
		break;
	case  CONTINUOUS_BMD_HYBRID_EXTRA:
	default:
		return -1;  // NO PARAMETER IS REMOVED THUS IT IS A NEGATIVE
					// INDEX AND WILL THROW AN ERROR
		break;
	}
}


/*****************************
* Function: bmd_start_extra_absolute
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the 
			closest point to the supplied value (usually the previous MLE)
			that satisfies the equality constraint. Note it always modifies parameters that are essentially
			unbounded - nu
*/
double  lognormalHILL_BMD_NC::bmd_start_absolute(unsigned n,
											  const double *b,
											  double *grad,
											  void   *data){
	
	start_data 		 *sdata = (start_data *) data; 
	Eigen::MatrixXd  theta  = sdata->theta; 
	
	if (!sdata->isIncreasing)
		sdata->BMRF *= -1.0; 
	
	double returnV = 0.0; 
	
	double gamma = theta(0,0); 
	double nu 	 = theta(1,0); 
	double k  	 = theta(2,0); 
	double n_exp = theta(3,0); 
	double temp  = 0.0; 
	returnV += pow(gamma - b[0], 2);
	returnV 	+= pow(k - b[2],2); 
	returnV 	+= pow(n_exp - b[3],2); 
	returnV		+= pow(theta(4, 0) - b[4], 2);
	temp 		 = sdata->BMRF*(pow(b[2],b[3])+pow(sdata->BMD,b[3])); 
	temp 		 = temp/pow(sdata->BMD,b[3]); 
	returnV 	+= pow(temp - nu,2.0); 

	return  returnV; 
}

std::vector<double> lognormalHILL_BMD_NC::bmd_start_absolute_clean(std::vector<double> x,
													double BMRF, double BMD,
													bool isIncreasing) 
{
	if (!isIncreasing)
		BMRF *= -1.0;
	double temp = BMRF; 
	temp = BMRF * (pow(x[2], x[3]) + pow(BMD, x[3]));
	temp = temp / pow(BMD, x[3]); 
	x[1] = temp; 
	return x; 
}

/**
* Function: bmd_start_reldev
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  lognormalHILL_BMD_NC::bmd_start_reldev(unsigned n,
											const double *b,
											double *grad,
											void   *data) {

	
	start_data 		 *sdata = (start_data *)data;
	Eigen::MatrixXd  theta = sdata->theta;

	Eigen::MatrixXd d(2,1); d << 0.0, sdata->BMD; 
	Eigen::MatrixXd theta_2 = theta; 
	for (int i = 0; i < n; i++) {theta_2(i,0) = b[i];}
	Eigen::MatrixXd t_median = mean(theta_2,d);
	t_median = exp(t_median.array()); 

	double t;

	if (sdata->isIncreasing) 
		t = sdata->BMRF;
	else
		t = 1.0 - sdata->BMRF;
	
	double returnV = 0.0;

	double gamma = theta(0, 0);
	double nu = theta(1, 0);
	double k = theta(2, 0);
	double n_exp = theta(3, 0);
	double temp = 0.0;
	
	returnV += pow(nu - b[1], 2);
	returnV += pow(k - b[2], 2);
	returnV += pow(n_exp - b[3], 2);
	returnV += pow(theta(4, 0) - b[4], 2);
	
	temp = fabs(t_median(1,0) - t_median(0,0))/t; 
	
	returnV += pow(temp - gamma, 2.0);


	return  returnV;
}

std::vector<double> lognormalHILL_BMD_NC::bmd_start_reldev_clean(std::vector<double> x,
	double BMRF, double BMD,
	bool isIncreasing)
{
	double t;
	Eigen::MatrixXd d(2,1); d << 0.0, BMD; 
	Eigen::MatrixXd theta_2(x.size(),1); ; 
	for (int i = 0; i < x.size(); i++) {theta_2(i,0) = x[i];}
	Eigen::MatrixXd t_median = mean(theta_2,d);
	t_median = exp(t_median.array());  

	if (isIncreasing)
		t = BMRF;
	else
		t = 1.0 - BMRF;

	double temp; 
	temp = fabs(x[1] * pow(BMD,x[3]) / (pow(x[2], x[3])+ pow(BMD, x[3])))/t; 
	x[0] = temp;
	return x;
}


/**********************************************
* Function: bmd_start_extra_relative
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
************************************************/
double  lognormalHILL_BMD_NC::bmd_start_stddev(unsigned n,
											const double *b,
											double *grad,
											void   *data) {
	// worry about the standard deviation here 
	start_data 		 *sdata = (start_data *)data;
	
	if (!sdata->isIncreasing)
		sdata->BMRF *= -1;
		
	Eigen::MatrixXd  theta = sdata->theta;
	Eigen::MatrixXd  theta_2 = theta; 
	
	for (int i = 0; i < n; i++) theta_2(i,0) = b[i]; 
	
	Eigen::MatrixXd d(2,1); d << 0.0, sdata->BMD; 
	Eigen::MatrixXd mu = mean(theta_2,d); 
	mu = exp(mu.array()); 
	double temp  = mu(1,0)/mu(0,0); 
	temp = 2.0*log(log(temp)/sdata->BMRF);  
	
	// squared Euclidean distance for the remaining parameters 
	double returnV = pow(temp-theta(n-1,0),2.0);  
	for (int i = 0; i < n-1; i++) returnV += pow( b[i] - theta(i,0),2.0); 
	
	return  returnV; 

}

std::vector<double> lognormalHILL_BMD_NC::bmd_start_stddev_clean(std::vector<double> x,
															  double BMRF, double BMD,
															  bool isIncreasing)
{
	if (!isIncreasing)
		BMRF *= -1;
	
	Eigen::MatrixXd  theta_2(x.size(),1); 
	
	for (int i = 0; i < x.size(); i++) theta_2(i,0) = x[i]; 
	
	Eigen::MatrixXd d(2,1); d << 0.0, BMD; 
	Eigen::MatrixXd mu = mean(theta_2,d); 
	
	double temp = fabs(x[1]*pow(BMD,x[3])/(pow(BMD,x[3])+pow(x[2],x[3])))/ x[0] + 1;
	temp = log(temp) / BMRF;
	temp = 2.0*log(temp);
	 
	x[x.size() - 1] = temp; 
	return x;
}

/****************************************************
* Function: bmd_start_extra_relative
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
***********************/
double  lognormalHILL_BMD_NC::bmd_start_point(unsigned n,
											const double *b,
											double *grad,
											void   *data) {

	start_data 		 *sdata = (start_data *)data;
	Eigen::MatrixXd  theta = sdata->theta;

	double returnV = 0.0;

	double gamma = theta(0, 0);
	double nu = theta(1, 0);
	double k = theta(2, 0);
	double n_exp = theta(3, 0);
	double temp = 0.0;

	returnV += pow(nu - b[1], 2);
	returnV += pow(k - b[2], 2);
	returnV += pow(n_exp - b[3], 2);
	temp += sdata->BMRF - (b[1] * pow(sdata->BMD , b[3])) / (pow(b[2],b[3]) + pow(b[2],b[3])); 
	returnV += pow(temp - gamma, 2.0); 
	
	temp = (pow(b[2], b[3]) + pow(sdata->BMD, b[3]));
	temp = b[1] * pow(sdata->BMD, b[3]) / temp * (1.0 / sdata->BMRF);
	temp = pow(temp, 2.0);
	temp = log(temp);
	returnV += pow(temp - b[4], 2.0);


	return  returnV;
}

std::vector<double> lognormalHILL_BMD_NC::bmd_start_point_clean(std::vector<double> x,
	double BMRF, double BMD,
	bool isIncreasing)
{
	double temp = 0;
	temp += BMRF - (x[1] * pow(BMD, x[3])) / (pow(x[2], x[3]) 
											+ pow(x[2], x[3]));
	x[0] = temp; 
	return x;
}


/**********************************************************
/*Function: bmd_start_extra
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*****************************/
double  lognormalHILL_BMD_NC::bmd_start_extra(unsigned n,
											const double *b,
											double *grad,
											void   *data) {

	start_data 		 *sdata = (start_data *)data;
	Eigen::MatrixXd  theta = sdata->theta;

	

	double returnV = 0.0;

	double gamma = theta(0, 0); 	double nu = theta(1, 0);
	double k = theta(2, 0);     	double n_exp = theta(3, 0);
	
		
	returnV += pow(k - b[2], 2);
	returnV += pow(n_exp - b[3], 2);
	returnV += pow(theta(4, 0) - b[4], 2);
	returnV += pow(b[1] - nu, 2.0);

	double temp = 0.0;
		temp = -1.0 / (sdata->BMRF)*b[1] * (pow(sdata->BMD, b[3]));
		temp = b[1] + temp / (pow(b[2], b[3]) + pow(sdata->BMD, b[3]));
	
	returnV += pow(gamma - temp, 2);

	return  returnV;
}

std::vector<double> lognormalHILL_BMD_NC::bmd_start_extra_clean(std::vector<double> x,
																double BMRF, double BMD,
																bool isIncreasing)
{
	
	double temp = 0.0;
	temp = -1.0 / (BMRF)*x[1] * (pow(BMD, x[3]));
	temp = x[1] + temp / (pow(x[2], x[3]) + pow(BMD, x[3]));
	

	x[0] = temp;
	return x;
}


/*********************************************************************
* Function: bmd_start_hybrid_extra
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  lognormalHILL_BMD_NC::bmd_start_hybrid_extra(unsigned n,
												  const double *b,
												  double *grad,
												  void   *data) {

	start_data 	*sdata = (start_data *)data;
	double	NOT_ADVERSE_P = 1.0 - sdata->tail_prob;
	double  TAIL_PROB  = sdata->tail_prob; 
	Eigen::MatrixXd  theta = sdata->theta;
	Eigen::MatrixXd	 theta2 = theta;
	/////////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < n; i++) { theta2(i, 0) = b[i]; }
	/////////////////////////////////////////////////////////////////////////////////////
	Eigen::MatrixXd d(2, 1); d << 0, sdata->BMD; 
	Eigen::MatrixXd temp_mean = mean(theta2, d);
	Eigen::MatrixXd temp_var = variance(theta2, d);
	/////////////////////////////////////////////////////////////////////////////////////
	double mu_zero = temp_mean(0, 0); double std_zero = sqrt(temp_var(0, 0));
	double ct_off = gsl_cdf_lognormal_Pinv(sdata->isIncreasing ? NOT_ADVERSE_P : 
																TAIL_PROB, mu_zero, std_zero);
	/////////////////////////////////////////////////////////////////////////////////////
	double returnV = 0.0; 
	double temp; 
	double k1 = gsl_cdf_ugaussian_Pinv(NOT_ADVERSE_P*sdata->BMRF+TAIL_PROB); 
	double k0 = gsl_cdf_ugaussian_Pinv(TAIL_PROB); 
	
    // Need to differenciate between increasing and decreasing
	if (sdata->isIncreasing){

		temp = (temp_mean(1,0) - temp_mean(0,0))/(k1-k0); // This is the standard deviation  
		temp = 2.0*log(temp); // transform it to the log scale and make it a variance
			
	}else{
		temp = (temp_mean(1, 0) - temp_mean(0, 0))/(k0-k1); // This is the standard deviation
		temp = 2.0*log(temp); // transform it to the log scale and make it a variance	
	}
	
	//////////////////////////////////////////////////////////////////////
	for (int i = 0; i <= n - 2; i++) {
		returnV += pow(theta(i,0) - b[i], 2.0);
	}
	returnV += pow(temp - theta(n-1,0), 2.0);	
	
	return  returnV;
	
}

std::vector<double> lognormalHILL_BMD_NC::bmd_start_hybrid_extra_clean(std::vector<double> x,
																double BMRF, double BMD,
																bool isIncreasing, double tail_prob)
{
	double	NOT_ADVERSE_P = 1.0 - tail_prob;
	double  TAIL_PROB = tail_prob;
	Eigen::MatrixXd	 theta2(x.size(),1);
	/////////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < x.size(); i++) { theta2(i, 0) = x[i]; }
	/////////////////////////////////////////////////////////////////////////////////////
	Eigen::MatrixXd d(2, 1); d << 0, BMD;
	Eigen::MatrixXd temp_mean = mean(theta2, d);
	Eigen::MatrixXd temp_var = variance(theta2, d);
	/////////////////////////////////////////////////////////////////////////////////////
	double mu_zero = temp_mean(0, 0); double std_zero = sqrt(temp_var(0, 0));
	double ct_off = gsl_cdf_lognormal_Pinv(isIncreasing ? NOT_ADVERSE_P :
											TAIL_PROB, mu_zero, std_zero);
	/////////////////////////////////////////////////////////////////////////////////////
	double returnV = 0.0;
	double temp;
	double k1 = gsl_cdf_ugaussian_Pinv(NOT_ADVERSE_P*BMRF + TAIL_PROB);

	double k0 = gsl_cdf_ugaussian_Pinv(TAIL_PROB);

	// Need to differenciate between increasing and decreasing
	if (isIncreasing) {

		temp = (temp_mean(1, 0) - temp_mean(0, 0)) / (k1 - k0); // This is the standard deviation  
		//cout << temp << endl; 
		temp = 2.0*log(temp); // transform it to the log scale and make it a variance

	}
	else {
		temp = (temp_mean(1, 0) - temp_mean(0, 0)) / (k0 - k1); // This is the standard deviation
		temp = 2.0*log(temp); // transform it to the log scale and make it a variance	
	}

	//////////////////////////////////////////////////////////////////////
	x[x.size()-1] = temp; // last term is ALWAYS the variance
	return x; 
}



/*************************************************
// Functions: lognormalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  lognormalHILL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  lognormalHILL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//            lognormalHILL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing) 
// Purpose :  return the BMD given the parameter values theta and the BMRF. Note they are  call the code
//            in  lognormalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
// 
**************************************************/
double lognormalHILL_BMD_NC::bmd_absolute_bound(Eigen::MatrixXd theta, double BMD,double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(2,1); d << 0.0, BMD; 
	Eigen::MatrixXd temp = mean(theta,d); 
	temp = exp(temp.array()); 
	
	double rValue = fabs(temp(0,0) - temp(1,0)) - BMRF; 
	return rValue; 
	
}

double lognormalHILL_BMD_NC::bmd_stdev_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = variance(theta,d); 
	Eigen::MatrixXd m	 =	  mean(theta, d);
	double t = temp(0,0); 
	
	return bmd_absolute_bound(theta,BMD,fabs(exp(m(0,0) + BMRF*pow(t,0.5))-exp(m(0,0))),isIncreasing); 	
	
}	

double lognormalHILL_BMD_NC::bmd_reldev_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	temp = exp(temp.array()); 
	
	double t; 
	
	if (isIncreasing)
		t = BMRF*temp(0,0); 
	else
		t = temp(0,0) - BMRF*temp(0,0); 
		
	
	return bmd_absolute_bound(theta,BMD,t,isIncreasing);
}

double lognormalHILL_BMD_NC::bmd_extra_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	temp = exp(temp.array()); 
	double nu = theta(1,0); // background mean
		
	if(isIncreasing)
		return bmd_absolute_bound(theta,BMD,BMRF*(nu-temp(0,0)),isIncreasing);   
	else
		return bmd_absolute_bound(theta,BMD,BMRF*(temp(0,0)-nu),isIncreasing);  
	
	
}

double lognormalHILL_BMD_NC::bmd_point_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	
	// note isIncreasing is ignored as point is specified by the user
	Eigen::MatrixXd d(1,1); d(0,0) = BMD; 
	Eigen::MatrixXd temp = mean(theta,d); 
	temp = exp(temp.array()); 
	
	return temp(0,0) - BMRF; 
}


/************************************************
////////////////////////////////////////////////////////////////////////
// Function:  double lognormalHILL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,,double BPROB)
// Variables: theta - matrix of theta values for the model
//			  BMRF  - This is a value between 0 and 1 that describes the increased probability over BPROB
//            isIncreasing - is the function an Increasing function or decreasing function? 
//            BPROB - Background probability at dose 0 considered adverse
// Purpose:   Compute the Hybrid BMD version of the hill model
//
//
//
////////////////////////////////////////////////////////////////////////
***********************************************/
double lognormalHILL_BMD_NC::bmd_hybrid_extra_bound(Eigen::MatrixXd theta, double BMD, 
										   double BMRF, bool isIncreasing,
										   double TAIL_PROB){
											   
											   
	///////////////////////////////////////////////////////////////////////////////////////
	double	NOT_ADVERSE_P = 1.0 - TAIL_PROB;
	Eigen::MatrixXd d(2, 1); d << 0.0, BMD;
	Eigen::MatrixXd mu  = mean(theta, d); // compute the mean at background an BMD
	Eigen::MatrixXd var = variance(theta, d);

	///////////////////////////////////////////////////////////////////////////////////////
	double ct_off = gsl_cdf_lognormal_Pinv(isIncreasing ? NOT_ADVERSE_P : TAIL_PROB, mu(0,0), sqrt(var(0,0)));
	///////////////////////////////////////////////////////////////////////////////////////
	//cout << ct_off << endl; 
	double temp; 
	
	if (isIncreasing){
		  temp = ((1.0-gsl_cdf_lognormal_P(ct_off,mu(1,0),sqrt(var(1,0))))-TAIL_PROB)/ NOT_ADVERSE_P;
 	}else{
		  temp = (gsl_cdf_lognormal_P(ct_off, mu(1, 0), sqrt(var(1, 0)))-TAIL_PROB)/ NOT_ADVERSE_P;
 	}

	return log(temp) - log(BMRF); 
}

/*********************************************************************************
// Functions: lognormalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  lognormalHILL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  lognormalHILL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//            lognormalHILL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing) 
// Purpose :  return the BMD given the parameter values theta and the BMRF. Note they are  call the code
//            in  lognormalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
// 
***********************************************************************************/
double lognormalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	
	if (!isIncreasing)
		BMRF *= -1.0; 
	
	double gamma = theta(0,0); 
	double nu = theta(1,0); 
	double k  = theta(2,0); 
	double n_exp = theta(3,0); 
	double temp = 0.0; 

	temp = (pow(BMRF/nu,1.0/n_exp))*k; 
	temp = temp/pow(1.0-BMRF/nu,1.0/n_exp); 

	return temp; 
	
}

double lognormalHILL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = variance(theta,d); 
	Eigen::MatrixXd med =  mean(theta, d);
	med = exp(med.array()); 
	double t = temp(0,0); 
	Eigen::MatrixXd md = abs(exp(log(med.array()) + BMRF * pow(t, 0.5)) - med.array()); 
	return bmd_absolute(theta,md(0,0),isIncreasing); 	
	
}	 

double lognormalHILL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	temp = exp(temp.array()); 
	
	double t; 
	
	if (isIncreasing)
		t = BMRF*temp(0,0); 
	else
		t = temp(0,0)*(1.0 - BMRF); 
		
	 
	return bmd_absolute(theta,t,isIncreasing);
}

double lognormalHILL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	// note isIncreasing is ignored as point is specified by the user
	
	double gamma = theta(0,0); 
	double nu = theta(1,0); 
	double k  = theta(2,0); 
	double n_exp = theta(3,0); 
	double temp = BMRF - gamma; 
	temp = temp/nu; 
	temp = temp/(1-temp); 
	
	return pow(temp,1.0/n_exp)*k; 
}

double lognormalHILL_BMD_NC::bmd_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	double nu = theta(1,0);
	
	if(isIncreasing)
		return bmd_absolute(theta,BMRF*(nu-temp(0,0)),isIncreasing);   
	else
		return bmd_absolute(theta,BMRF*(temp(0,0)-nu),isIncreasing);   
}

/**********************************************************************
////////////////////////////////////////////////////////////////////////
// Function:  double lognormalHILL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,,double BPROB)
// Variables: theta - matrix of theta values for the model
//			  BMRF  - This is a value between 0 and 1 that describes the increased probability over BPROB
//            isIncreasing - is the function an Increasing function or decreasing function? 
//            BPROB - Background probability at dose 0 considered adverse
// Purpose:   Compute the Hybrid BMD version of the hill model
//
//
//
////////////////////////////////////////////////////////////////////////
**************************************************************************/
double lognormalHILL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,double TAIL_PROB){
    
	double NOT_ADVERSE_P = 1.0 - TAIL_PROB;

	////////////////////////////////////////////////////////////////////
	//Get the mean and variance at dose zero as well as a very high dose
	double min_d = 0.0; double max_d =  X.maxCoeff(); double mid = 0.5*(min_d + max_d);
	Eigen::MatrixXd d(3, 1); d << min_d, mid, max_d;
	Eigen::MatrixXd temp_mean = mean(theta, d);
	Eigen::MatrixXd temp_var = variance(theta, d);
	//////////////////////////////////////////////////////////////////////
	double mu_zero = temp_mean(0, 0); double std_zero = sqrt(temp_var(0, 0));
	double ct_off = gsl_cdf_lognormal_Pinv(isIncreasing ? NOT_ADVERSE_P : TAIL_PROB, mu_zero, std_zero); // CUTOFF AT DOSE = 0
	double P = TAIL_PROB + BMRF * NOT_ADVERSE_P;
	//double bmr_mult =  0.0; //gsl_cdf_lognormal_Pinv(NOT_ADVERSE_P - BMRF * (NOT_ADVERSE_P));

	double test_prob = isIncreasing ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)))
		: gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0))); //standardize
	double test = 0;


	int k = 0;
	while (test_prob < P  && k < 10) { // Go up to 2^10 times the maximum tested dose
									   // if we cant find it after that we return infinity
		max_d *= 2;
		d << min_d, mid, max_d;
		temp_mean = mean(theta, d);
		temp_var = variance(theta, d);

		test_prob = isIncreasing ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)))
			: gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)));

		k++;
	}

	if (k == 10 || test_prob < P) // have not been able to bound the BMD
	{
		return std::numeric_limits<double>::infinity();
	}


	test_prob = isIncreasing ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)))
		: gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)));

	double temp_test = test_prob - P;
	///////////////////////////////////////////////////////////////////////////// 
	while (fabs(temp_test) > 1e-5) {
		// we have bounded the BMD now we use a root finding algorithm to 
		// figure out what it is default difference is a probability of of 1e-5
		if (temp_test > 0) {
			max_d = mid;
		}
		else {
			min_d = mid;
		}

		mid = 0.5*(max_d + min_d);
		d << min_d, mid, max_d;

		temp_mean = mean(theta, d);
		temp_var = variance(theta, d);

		test_prob = isIncreasing ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)))
			: gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)));

		temp_test = test_prob - P;
	}


	if (isfinite(mid)){
	  return mid; 
	}else{
	 // std::cerr << "Non-finite BMD returned: Exp-Log-Hill."<< std::endl;
	  return std::numeric_limits<double>::infinity();
	}
}
