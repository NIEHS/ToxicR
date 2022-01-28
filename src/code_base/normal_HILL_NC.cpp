#include "stdafx.h" // Precompiled header - does nothing if building R version
#include "normal_HILL_NC.h"
#include "normalModels.h"

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
    
#endif

#include <gsl/gsl_randist.h>
#include <iostream>
#include <cmath>
#include <math.h>  
using namespace std;

/////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////
int normalHILL_BMD_NC::parameter_to_remove(contbmd TYPE) {

	switch (TYPE) {
	case CONTINUOUS_BMD_ABSOLUTE:
		return 1; // slope
	  break; 
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


/*Function: bmd_start_extra_absolute
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the 
			closest point to the supplied value (usually the previous MLE)
			that satisfies the equality constraint. Note it always modifies parameters that are essentially
			unbounded - nu
*/
double  normalHILL_BMD_NC::bmd_start_absolute(unsigned n,
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
	if (n == 6) {
		returnV += pow(theta(5, 0) - b[5], 2.0); 
	}
	
	return  returnV; 
}

std::vector<double> normalHILL_BMD_NC::bmd_start_absolute_clean(std::vector<double> x,
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

/*Function: bmd_start_reldev
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  normalHILL_BMD_NC::bmd_start_reldev(unsigned n,
											const double *b,
											double *grad,
											void   *data) {

	
	start_data 		 *sdata = (start_data *)data;
	Eigen::MatrixXd  theta = sdata->theta;

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
	temp = (pow(b[2], b[3]) + pow(sdata->BMD, b[3]));
	temp = b[1]*pow(sdata->BMD, b[3])/temp * (1.0/t);
	temp = sdata->isIncreasing ? temp : -1.0*temp;
	returnV += pow(temp - gamma, 2.0);
	
	if (n == 6) {
		returnV += pow(theta(5, 0) - b[5], 2.0);
	}

	return  returnV;
}

std::vector<double> normalHILL_BMD_NC::bmd_start_reldev_clean(std::vector<double> x,
	double BMRF, double BMD,
	bool isIncreasing)
{
	double t;

	if (isIncreasing)
		t = BMRF;
	else
		t = 1.0 - BMRF;

	double temp; 
	temp = (pow(x[2], x[3]) + pow(BMD, x[3]));
	temp = x[1] * pow(BMD, x[3]) / temp * (1.0 / t);
	
	x[0] = isIncreasing? temp:-1.0*temp;
	return x;
}


/*Function: bmd_start_extra_relative
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  normalHILL_BMD_NC::bmd_start_stddev(unsigned n,
											const double *b,
											double *grad,
											void   *data) {


	start_data 		 *sdata = (start_data *)data;
	if (!sdata->isIncreasing)
		sdata->BMRF *= 1;

	Eigen::MatrixXd  theta = sdata->theta;
	Eigen::MatrixXd  theta_2 = theta;

	for (int i = 0; i < n; i++) theta_2(i, 0) = b[i];

	Eigen::MatrixXd d(2, 1); d << 0.0, sdata->BMD;
	Eigen::MatrixXd mu = mean(theta_2, d);

	double temp;

	if (!isConstVar()) {
		temp = log(fabs((mu(1, 0) - mu(0, 0))));
		temp -= log(sdata->BMRF);
		temp -= log(mu(0, 0))*b[n - 2] / 2.0;
	}
	else {
		temp = log(fabs((mu(1, 0) - mu(0, 0))));
		temp -= log(sdata->BMRF);
	}
	temp = 2.0*temp;

	double returnV = pow(temp - theta(n - 1, 0), 2.0);

	// squared Euclidean distance for the remaining parameters 
	for (int i = 0; i < n - 1; i++) returnV += pow(b[i] - theta(i, 0), 2.0);

	return  returnV;

}

std::vector<double> normalHILL_BMD_NC::bmd_start_stddev_clean(std::vector<double> x,
															  double BMRF, double BMD,
															  bool isIncreasing)
{ 
	if (!isIncreasing)
		BMRF *= 1;


	Eigen::MatrixXd  theta_2(x.size(), 1);
	//cout << BMD << " "; 
	for (int i = 0; i < x.size(); i++)  theta_2(i, 0) = x[i];
	//for (double b : x) cout << b << " "; 
	Eigen::MatrixXd d(2, 1); d << 0.0, BMD;
	Eigen::MatrixXd mu = mean(theta_2, d);

	double temp;

	if (!isConstVar()) {
		temp = log(fabs((mu(1, 0) - mu(0, 0))));
		temp -= log(BMRF) + log(mu(0, 0))*x[x.size() - 2] / 2.0;
	}
	else {
		temp = log(fabs((mu(1, 0) - mu(0, 0))));
		temp -= log(BMRF);
	}
	temp = 2.0*temp;
 
	x[x.size() - 1] = temp;
	return x;
}

/*Function: bmd_start_extra_relative
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  normalHILL_BMD_NC::bmd_start_point(unsigned n,
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
	
	if (n == 6) { // non-constant variance
		returnV += pow(theta(4, 0) - b[4], 2.0);
		returnV += pow(theta(5, 0) - b[5], 2.0);
	}
	else {  // constant variance
		returnV += pow(theta(4, 0) - b[4], 2.0);
	}

	return  returnV;
}

std::vector<double> normalHILL_BMD_NC::bmd_start_point_clean(std::vector<double> x,
	double BMRF, double BMD,
	bool isIncreasing)
{
	double temp = 0;
	temp += BMRF - (x[1] * pow(BMD, x[3])) / (pow(x[2], x[3]) 
											+ pow(x[2], x[3]));
	x[0] = temp; 
	return x;
}



/*Function: bmd_start_extra
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  normalHILL_BMD_NC::bmd_start_extra(unsigned n,
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
	if (n == 6) {
		returnV += pow(theta(5, 0) - b[5], 2.0);
	}

	return  returnV;
}

std::vector<double> normalHILL_BMD_NC::bmd_start_extra_clean(std::vector<double> x,
																double BMRF, double BMD,
																bool isIncreasing)
{
	
	double temp = 0.0;
	temp = -1.0 / (BMRF)*x[1] * (pow(BMD, x[3]));
	temp = x[1] + temp / (pow(x[2], x[3]) + pow(BMD, x[3]));
	

	x[0] = temp;
	return x;
}


/*Function: bmd_start_hybrid_extra
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double  normalHILL_BMD_NC::bmd_start_hybrid_extra(unsigned n,
												  const double *b,
												  double *grad,
												  void   *data) {

	start_data 	*sdata = (start_data *)data;
	double	BPROB = 1.0 - sdata->tail_prob; 
	double  TAIL  = sdata->tail_prob; 
	Eigen::MatrixXd  theta = sdata->theta;



	double returnV = 0.0;

	double gamma = theta(0, 0); 	double nu = theta(1, 0);
	double k = theta(2, 0);     	double n_exp = theta(3, 0);
	
		
	

	Eigen::MatrixXd d(2,1); d << 0.0, sdata->BMD; 
	Eigen::MatrixXd thetas_d(n,1); 
	for (int i = 0; i < n; i++){ thetas_d(i,0) = b[i];};
	Eigen::MatrixXd mu = mean(thetas_d,d); // compute the mean at background an BMD
	

	//////////////////////////////////////////////////////////////////////
    
    double k_1 = gsl_cdf_ugaussian_Pinv(BPROB*sdata->BMRF+TAIL); 
    double k_0 = gsl_cdf_ugaussian_Pinv(TAIL); 
    double l; 
    double temp; 
    
    // Need to differenciate between increasing and decreasing
	if (sdata->isIncreasing){
		l = mu(1,0) - mu(0,0);
		if (n == 6)
			temp = pow(mu(1,0),b[4]/2.0)*k_1 - pow(mu(0,0),b[4]/2.0)*k_0 ;  
		else
			temp = k_1 - k_0; 
		
		l = l/temp;		   // This is the standard deviation
		temp = 2.0*log(l); // transform it to the log scale and make it a variance
		
	}else{
		
		l = mu(1,0) - mu(0,0);
		if (n == 6)
			temp = pow(mu(0,0),b[4]/2.0)*k_0  - pow(mu(1,0),b[4]/2.0)*k_1 ;  
		else
			temp = k_1 - k_0; 
		
		l = l/temp;		   // This is the standard deviation
		temp = 2.0*log(l); // transform it to the log scale and make it a variance	
	}
	//////////////////////////////////////////////////////////////////////


	returnV += pow(gamma - b[0], 2.0);
	returnV += pow(b[1] - nu, 2.0);
	returnV += pow(k - b[2], 2.0);
	returnV += pow(n_exp - b[3], 2.0);
	
	if (n == 6){ // power model
		returnV += pow(theta(4, 0) - b[4], 2.0); 
		returnV += pow(theta(5, 0) - temp, 2.0); 
	}else{		
		returnV += pow(temp - theta(4,0), 2.0);	
	}
	

	return  returnV;
}

std::vector<double> normalHILL_BMD_NC::bmd_start_hybrid_extra_clean(std::vector<double> x,
																double BMRF, double BMD,
																bool isIncreasing, double tail_prob)
{
	
	
	double	BPROB = 1.0 - tail_prob; 
	double  TAIL  = tail_prob; 
	
	

	Eigen::MatrixXd d(2,1); d << 0.0, BMD; 
	Eigen::MatrixXd thetas_d(x.size(),1); 
	for (int i = 0; i < x.size(); i++){ thetas_d(i,0) = x[i];};
	Eigen::MatrixXd mu = mean(thetas_d,d); // compute the mean at background an BMD
	

	//////////////////////////////////////////////////////////////////////
    
    double k_1 = gsl_cdf_ugaussian_Pinv(BPROB*BMRF+TAIL); 
    double k_0 = gsl_cdf_ugaussian_Pinv(TAIL); 
    double l; 
    double temp; 
    
    // Need to differenciate between increasing and decreasing
	if (isIncreasing){
		l = mu(1,0) - mu(0,0);
		if (x.size() == 6)
			temp = pow(mu(1,0),x[4]/2.0)*k_1 - pow(mu(0,0),x[4]/2.0)*k_0  ; 
		else
			temp = k_1 - k_0; 
		
		l = l/temp;		   // This is the standard deviation
		temp = 2.0*log(l); // transform it to the log scale and make it a variance
		
	}else{
		
		l = mu(1,0) - mu(0,0);
		if (x.size() == 6)
			temp = pow(mu(0,0),x[4]/2.0)*k_0  - pow(mu(1,0),x[4]/2.0)*k_1  ; 
		else
			temp = k_0 - k_1; 
		
		l = l/temp;		   // This is the standard deviation
		temp = 2.0*log(l); // transform it to the log scale and make it a variance	
	}
	//////////////////////////////////////////////////////////////////////
	
	if (x.size() == 6){ // power model
		x[5] = temp; 
	}else{		
		x[4] = temp; 
	}

	return x; 
}




// Functions: normalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  normalHILL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  normalHILL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//            normalHILL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing) 
// Purpose :  return the BMD given the parameter values theta and the BMRF. Note they are  call the code
//            in  normalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
// 
double normalHILL_BMD_NC::bmd_absolute_bound(Eigen::MatrixXd theta, double BMD,double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(2,1); d << 0.0, BMD; 
	Eigen::MatrixXd temp = mean(theta,d); 
	
	
	double rValue = fabs(temp(0,0) - temp(1,0)) - BMRF; 
	return rValue; 
	
}

double normalHILL_BMD_NC::bmd_stdev_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = variance(theta,d); 
	double t = temp(0,0); 
	
	return bmd_absolute_bound(theta,BMD,BMRF*pow(t,0.5),isIncreasing); 	
	
}	 

double normalHILL_BMD_NC::bmd_reldev_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	
	double t; 
	
	if (isIncreasing)
		t = BMRF*temp(0,0); 
	else
		t = temp(0,0) - BMRF*temp(0,0); 
		
	
	return bmd_absolute_bound(theta,BMD,t,isIncreasing);
}

double normalHILL_BMD_NC::bmd_extra_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	double nu = theta(1,0); // background mean
		
	if(isIncreasing)
		return bmd_absolute_bound(theta,BMD,BMRF*(nu-temp(0,0)),isIncreasing);   
	else
		return bmd_absolute_bound(theta,BMD,BMRF*(temp(0,0)-nu),isIncreasing);  
	
	
}

double normalHILL_BMD_NC::bmd_point_bound(Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing){
	
	// note isIncreasing is ignored as point is specified by the user
	Eigen::MatrixXd d(1,1); d(0,0) = BMD; 
	Eigen::MatrixXd temp = mean(theta,d); 
	
	return temp(0,0) - BMRF; 
}



////////////////////////////////////////////////////////////////////////
// Function:  double normalHILL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,,double BPROB)
// Variables: theta - matrix of theta values for the model
//			  BMRF  - This is a value between 0 and 1 that describes the increased probability over BPROB
//            isIncreasing - is the function an Increasing function or decreasing function? 
//            BPROB - Background probability at dose 0 considered adverse
// Purpose:   Compute the Hybrid BMD version of the hill model
//
//
//
////////////////////////////////////////////////////////////////////////
double normalHILL_BMD_NC::bmd_hybrid_extra_bound(Eigen::MatrixXd theta, double BMD, 
										   double BMRF, bool isIncreasing,
										   double TAIL_PROB){
											   
	////////////////////////////////////////////////////////////////////
	double	BPROB = 1.0 - TAIL_PROB; 
    ////////////////////////////////////////////////////////////////////
    //Get the mean and variance at dose zero as well as a very high dose
	double min_d = 0.0;
	Eigen::MatrixXd d(2,1); 		 d << min_d, BMD; 	
	Eigen::MatrixXd temp_mean 	=      mean(theta,d); 
	Eigen::MatrixXd temp_var    =  variance(theta,d); 
	
	double mu_zero = temp_mean(0,0); double std_zero = pow(temp_var(0,0),0.5); 
	double mu_bmd =  temp_mean(1,0); double std_bmd  = pow(temp_var(1,0),0.5); 
	
//	cout << mu_zero << " "<< std_zero << " " << 
//		    mu_bmd  << " "<<  std_bmd << " " << endl; 
	
	double l; 
	double temp; 
	
	if (isIncreasing){
		  l = mu_zero - gsl_cdf_ugaussian_Pinv(TAIL_PROB)*std_zero; 
		  temp = (gsl_cdf_gaussian_P( mu_bmd-l,std_bmd)  - TAIL_PROB)/BPROB; 
 	}else{
		  l = mu_zero + gsl_cdf_ugaussian_Pinv(TAIL_PROB)*std_zero; 
		  
		  temp = (gsl_cdf_gaussian_P( l-mu_bmd,std_bmd)  - TAIL_PROB)/BPROB; 
 	}
 	
 	//cout << temp << endl;  
	return log(temp) - log(BMRF); 
}

// Functions: normalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  normalHILL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//			  normalHILL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
//            normalHILL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing) 
// Purpose :  return the BMD given the parameter values theta and the BMRF. Note they are  call the code
//            in  normalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing)
// 
double normalHILL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	
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
double normalHILL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = variance(theta,d); 
	double t = temp(0,0); 
	
	return bmd_absolute(theta,BMRF*pow(t,0.5),isIncreasing); 	
	
}	 
double normalHILL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	
	double t; 
	
	if (isIncreasing)
		t = BMRF*temp(0,0); 
	else
		t = temp(0,0) - BMRF*temp(0,0); 
		
	 
	return bmd_absolute(theta,t,isIncreasing);
}

double normalHILL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
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

double normalHILL_BMD_NC::bmd_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing){
	Eigen::MatrixXd d(1,1); d(0,0) = 0.0; 
	Eigen::MatrixXd temp = mean(theta,d); 
	double nu = theta(1,0);
	
	if(isIncreasing)
		return bmd_absolute(theta,BMRF*(nu-temp(0,0)),isIncreasing);   
	else
		return bmd_absolute(theta,BMRF*(temp(0,0)-nu),isIncreasing);   
}

////////////////////////////////////////////////////////////////////////
// Function:  double normalHILL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,,double BPROB)
// Variables: theta - matrix of theta values for the model
//			  BMRF  - This is a value between 0 and 1 that describes the increased probability over BPROB
//            isIncreasing - is the function an Increasing function or decreasing function? 
//            BPROB - Background probability at dose 0 considered adverse
// Purpose:   Compute the Hybrid BMD version of the hill model
//
//
//
////////////////////////////////////////////////////////////////////////
double normalHILL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double BMRF, bool isIncreasing,double TAIL_PROB){
    
	
	double NOT_ADVERSE_P = 1.0 - TAIL_PROB;
    
  ////////////////////////////////////////////////////////////////////
  //Get the mean and variance at dose zero as well as a very high dose
	double min_d = 0.0; double max_d =  X.maxCoeff(); double mid = 0.5*(min_d+max_d); 
	Eigen::MatrixXd d(3,1); d << min_d, mid, max_d; 	
	Eigen::MatrixXd temp_mean =     mean(theta,d); 
	Eigen::MatrixXd temp_var  = variance(theta,d); 
	
	double mu_zero   = temp_mean(0,0); double std_zero = pow(temp_var(0,0),0.5); 
 	double inv_bprob = gsl_cdf_ugaussian_Pinv(NOT_ADVERSE_P);
 	
	double bmr_mult = gsl_cdf_ugaussian_Pinv(NOT_ADVERSE_P - BMRF * (NOT_ADVERSE_P));
									
		
 	double temp_const = (isIncreasing)?inv_bprob*std_zero+mu_zero: mu_zero-inv_bprob * std_zero ;
 	
    ////////////////////////////////////////////////////////////////////
	double test =  temp_const - temp_mean(2,0) ; //standardize
	test  *= 1/pow(temp_var(2,0),0.5); 
	double test_prob  = isIncreasing? gsl_cdf_ugaussian_P(test): 1.0 - gsl_cdf_ugaussian_P(test); 
	double P = gsl_cdf_ugaussian_P(bmr_mult); 
	
	int k = 0; 
	while (test_prob > P  && k <  10){ // Go up to 2^10 times the maximum tested dose
						         // if we cant find it after that we return infinity
	     max_d *= 2.0; 
	     mid = 0.5*(min_d+max_d); 
	     d << min_d, mid, max_d; 
	     temp_mean =     mean(theta,d); 
		 temp_var  = variance(theta,d); 
	
		test =  temp_const - temp_mean(2,0) ; //standardize
		test  *= 1/pow(temp_var(2,0),0.5); 
		
		test_prob = isIncreasing ?  gsl_cdf_ugaussian_P(test) :  1.0 - gsl_cdf_ugaussian_P(test);
		k ++; 
	}
	
 
	if (k == 10) // have not been able to bound the BMD
	{
		return std::numeric_limits<double>::infinity();
	}
	
	test =  temp_const - temp_mean(1,0) ; //standardize
	test  *= 1/pow(temp_var(1,0),0.5); 
	test_prob = isIncreasing ? gsl_cdf_ugaussian_P(test) : 1.0 - gsl_cdf_ugaussian_P(test);
	double temp_test = test_prob - P;
	 

  int iter = 0; 
	while (fabs(temp_test) > 1e-5 && iter < 200){
		// we have bounded the BMD now we use a root finding algorithm to 
		// figure out what it is default difference is a probability of of 1e-5
		if (temp_test  < 0){
				max_d = mid; 
		} else {
				min_d = mid; 
		}
		mid = 0.5*(max_d+min_d);
		d << min_d, mid, max_d; 
		
		temp_mean =     mean(theta,d); 
		temp_var  = variance(theta,d);  
			
		test =  temp_const - temp_mean(1,0) ; //standardize
		test  *= 1/pow(temp_var(1,0),0.5); 
		test_prob = isIncreasing ? gsl_cdf_ugaussian_P(test) : 1.0 - gsl_cdf_ugaussian_P(test);
					
		temp_test = test_prob - P;
    iter ++; 
		
	}
	
	if (isfinite(mid)){
	  return mid; 
	}else{
	  // std::cerr << "Non-finite BMD returned: Hill-Normal."<< std::endl; 
	  return std::numeric_limits<double>::infinity();
	}
}

