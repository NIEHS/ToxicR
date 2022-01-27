

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

    
#include <vector>
using namespace std; 
//////////////////////////////////////////////////////////////////
//function: standard_data
//purpose: takes data either origional or sufficient statistics
//makes everythign 
//
//
//////////////////////////////////////////////////////////////////
Eigen::MatrixXd standard_data(Eigen::MatrixXd Y,Eigen::MatrixXd X,double *zero_mean){
	
	Eigen::MatrixXd returnV = Y; 
	
	// Find the average of the zero dose group
	// This will be used to change the data.
	// This is for the Bayesian analysis, so the priors are consistant	
	// with %change from background. x
	vector <double> temp; 
	for (int i = 0; i < X.rows(); i++){
			if(X(i,0) == 0){
				temp.push_back(Y(i,0)); // find all the zero dose columns
			}
	}
	double mean = 0.0; 
	for (int i = 0; i < temp.size(); i++){
			mean += temp[i]/temp.size();  // find the mean of the zero 
	}
	
		
	for (int i = 0; i < returnV.rows(); i++){
			returnV(i,0) *= 1/mean; 
			if(returnV.cols()==3){ // sufficient statistics
					returnV(i,2) *= 1/mean; 
			}
	}
	*zero_mean = mean; 
	return returnV; 
}
