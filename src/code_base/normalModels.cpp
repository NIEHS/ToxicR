// File: normalModels.h
// Purpose: 
// Implements the normal model likelihoods that are used for BMD analysis
//
//  Creator:       Matt Wheeler
//  Creation Date: 4/17/2018
//
//
//

#include "stdafx.h" // Precompiled header - does nothing if building R version
#include "normalModels.h"

#ifdef R_COMPILATION
    //necessary things to run in R    starting_value
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif
#include "cmodeldefs.h"
#include "normal_likelihoods.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/////////////////////////////////////////////////////////////////////////
// function: type_of_profile()
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
int normalLLModel::type_of_profile(contbmd TYPE) {	
	
	switch (TYPE) {
		case CONTINUOUS_BMD_ABSOLUTE:
		case CONTINUOUS_BMD_STD_DEV:
		case CONTINUOUS_BMD_REL_DEV:
		case CONTINUOUS_BMD_POINT:
		case CONTINUOUS_BMD_EXTRA:
			return PROFILE_INEQUALITY;
			break;
		case  CONTINUOUS_BMD_HYBRID_EXTRA:
		default: 
			return PROFILE_EQUALITY;  // HYBRID IS ALWAYS AN EQUALITY CONSTRAINT
			break;
	}
}

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
int normalLLModel::parameter_to_remove(contbmd TYPE) {
	
	switch (TYPE) {
	case CONTINUOUS_BMD_ABSOLUTE:
	case CONTINUOUS_BMD_STD_DEV:
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

////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////

double start_distance(unsigned     n,
					 const double *b,
					 double    *grad,
					 void      *data) {
	start_data    *sdata = (start_data *)data;
	normalLLModel *model = (normalLLModel*) sdata->M;
	
	switch (sdata->BMDType) {
	case CONTINUOUS_BMD_ABSOLUTE:
		return model->bmd_start_absolute(n, b, grad, data); 
		break;
	case  CONTINUOUS_BMD_STD_DEV:
		return model->bmd_start_stddev(n, b, grad, data);
		break;
	case CONTINUOUS_BMD_REL_DEV:
		return model->bmd_start_reldev(n, b, grad, data);
		break;
	case CONTINUOUS_BMD_POINT:
		return model->bmd_start_point(n, b, grad, data); 
		break;
	case CONTINUOUS_BMD_EXTRA:
		return model->bmd_start_extra(n, b, grad, data);
		break;
	case  CONTINUOUS_BMD_HYBRID_EXTRA:
		return model->bmd_start_hybrid_extra(n, b, grad, data);
		break;
	case CONTINUOUS_BMD_HYBRID_ADDED:
		break;
	default:
		return 0.0; 
	}
	 return 0.0; 
	
}
///////////////////////////////////////////////////////////////////////
// Function: starting_value()
// Purpose : for equality constrained optimization finds a value of x 
//           that is on the equality constraint and satisfies the bounds. 
// Input: Origional value. 
// Output: Value that satisfies the equality constraint of the given
//         risk type.  
///////////////////////////////////////////////////////////////////////
Eigen::MatrixXd normalLLModel::starting_value(Eigen::MatrixXd theta, contbmd BMDType, double BMD, double BMRF,
												bool isInc,
												double tail_prob,
												std::vector<double> lb,
												std::vector<double> ub) 
{
	////////////////////////////////////////////
	double minf;
	nlopt::result result;
	////////////////////////////////////////////
	start_data data; 
	data.M = this; data.BMDType = BMDType; 
	data.theta = theta; data.BMDType = BMDType; data.BMD=BMD; 
	data.BMRF  = BMRF ; data.isIncreasing   = isInc; 
	data.tail_prob = tail_prob; 
	
	nlopt::opt opt(nlopt::LN_BOBYQA, theta.rows());
	opt.set_lower_bounds(lb); 
	opt.set_upper_bounds(ub);
	opt.set_xtol_abs(1e-5);
	opt.set_maxeval(20000);
	
	
	std::vector<double> x(theta.rows()); 
	for (int i = 0; i < theta.rows(); i++)
	   x[i] = theta(i, 0);
	
	std::vector<double> init(x.size());
	for (int i = 0; i < x.size(); i++) init[i] = 1e-4;
	opt.set_initial_step(init);
	
	opt.set_min_objective(start_distance,&data);
		
	bool good_opt = false; 
    result = nlopt::FAILURE;
	try{
			result = opt.optimize(x, minf);
			good_opt =  true;  
	}catch (nlopt::roundoff_limited &exec) {
			good_opt = false; 
	}
	catch (nlopt::forced_stop &exec) {
			good_opt = false;
	}
	catch (const std::exception &exc) {
			good_opt = false;	
	}
	if (result > 4){ // Either 5 = Max iterations or 6 = Max Time 
			good_opt = false; 
	}	
	
	switch (BMDType) {
	case CONTINUOUS_BMD_ABSOLUTE:
		x = bmd_start_absolute_clean(x, BMRF, BMD, isInc);
		break;
	case  CONTINUOUS_BMD_STD_DEV:
		x = bmd_start_stddev_clean(x, BMRF, BMD, isInc); 
		break;
	case CONTINUOUS_BMD_REL_DEV:
		x = bmd_start_reldev_clean(x, BMRF, BMD, isInc);
		break;
	case CONTINUOUS_BMD_POINT:
		x = bmd_start_point_clean(x, BMRF, BMD, isInc); 
		break;
	case CONTINUOUS_BMD_EXTRA:
		x = bmd_start_extra_clean(x, BMRF, BMD, isInc);
		break;
	case  CONTINUOUS_BMD_HYBRID_EXTRA:
		x = bmd_start_hybrid_extra_clean(x, BMRF, BMD, isInc,tail_prob);
		
		break;
	case CONTINUOUS_BMD_HYBRID_ADDED:
		break;
	default:
		break; 
	}

	
	Eigen::Map<Eigen::MatrixXd> d(x.data(), theta.rows(), 1);

	if (good_opt) {
		return d; 
	}else {
		return Eigen::MatrixXd::Zero(theta.rows(), 1);
	}

}
//////////////////////////////////////////////////////
//
//
//
//
//////////////////////////////////////////////////////////
double normalLLModel::equality_boundG(Eigen::MatrixXd theta, contbmd BMDType, double BMD, double BMRF,
												bool isInc,
												double tail_prob){
							   
		double returnV = 0.0; 	
		switch (BMDType){
			case CONTINUOUS_BMD_ABSOLUTE:    
				 returnV = bmd_absolute_bound( theta,BMD, BMRF,isInc); 
				 break; 
			case  CONTINUOUS_BMD_STD_DEV:   
			     returnV = bmd_stdev_bound( theta, BMD, BMRF,isInc);  
				break; 
			case CONTINUOUS_BMD_REL_DEV :    
				  returnV = bmd_reldev_bound( theta, BMD, BMRF,isInc); 
				break; 
			case CONTINUOUS_BMD_POINT   :   
				  returnV = bmd_point_bound( theta, BMD, BMRF,isInc);  
				break;  
			case CONTINUOUS_BMD_EXTRA   :    
				  returnV = bmd_extra_bound( theta, BMD, BMRF,isInc); 	
				break; 
			case  CONTINUOUS_BMD_HYBRID_EXTRA: 
				  returnV = bmd_hybrid_extra_bound(theta, BMD, BMRF,isInc,tail_prob); 
				break;
			case CONTINUOUS_BMD_HYBRID_ADDED:
				break; 
			default: 
				break; 	
		}

		return returnV; 
							   


}

///////////////////////////////////////////////////////////////////////
//Function: eqConst_gradient - Computes the gradient of the equality constraint
///////////////////////////////////////////////////////////////////////
Eigen::MatrixXd normalLLModel::eqConst_gradient(Eigen::MatrixXd v, contbmd type,
												double BMD, double BMRF,
												bool isIncreasing, double tail_prob) {

	Eigen::VectorXd h(v.rows());
	double mpres = pow(1.0e-16, 0.5);
	double x, temp;
	double derEst;
	double f1 = 0.0; double f2 = 0.0;
	Eigen::MatrixXd  hvector = v;
	Eigen::MatrixXd rValue(v.rows(), 1);

	for (int i = 0; i < v.rows(); i++)
	{
		x = v(i, 0);
		if (fabs(x) > DBL_EPSILON) {
			h[i] = mpres * (fabs(x));
			temp = x + h[i];
			h[i] = temp - x;
		}
		else {
			h[i] = mpres;
		}

	}

	/*find the gradients for each of the variables in the function*/
	for (int i = 0; i < v.rows(); i++) {
		/*perform a finite difference calculation on the specific derivative*/
		x = v(i, 0);
		// add h
		hvector(i, 0) = x + h[i];
		f1 = equality_boundG(hvector, type, BMD, BMRF, isIncreasing, tail_prob);
		// subtract h
		hvector(i, 0) = x - h[i];
		f2 = equality_boundG(hvector, type, BMD, BMRF, isIncreasing, tail_prob);
		// estimate the derivative
		rValue(i, 0) = (f1 - f2) / (2.0*h[i]);
		hvector(i, 0) = x;
	}

	return rValue;
}



Eigen::MatrixXd normalLLModel::mean(Eigen::MatrixXd theta) {
		return mean(theta,X); 
}										    							    	
	
Eigen::MatrixXd normalLLModel::mean(Eigen::MatrixXd theta, Eigen::MatrixXd d) {
		double beta0 = theta(0,0); 
		double beta1 = theta(1,0); 
		
		Eigen::MatrixXd rV = beta0 + beta1*d.array();		
		
		return rV; 
}

Eigen::MatrixXd normalLLModel::variance(Eigen::MatrixXd theta) {
		return variance(theta, X); 
}
	
Eigen::MatrixXd normalLLModel::variance(Eigen::MatrixXd theta, Eigen::MatrixXd d){
	  double var = exp(theta(theta.rows()-1,0)); 
		Eigen::MatrixXd rV = Eigen::MatrixXd::Ones(d.rows(),1)*var; 
		 
		if (!isConstVar()){ //variance is proportional to power function of the mean
			Eigen::MatrixXd zero_dose(1, 1); zero_dose << 0;
			Eigen::MatrixXd mean_zero = mean(theta, d); 
			Eigen::MatrixXd meanV = mean(theta,d);
			rV = exp(log( abs(meanV.array()))*theta(theta.rows()-2,0));		
			rV = var * rV;
		}
		 
		return rV;  					
}
