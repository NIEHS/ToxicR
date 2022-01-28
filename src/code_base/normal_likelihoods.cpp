#include <normal_likelihoods.h>


double normalLL::negLogLikelihood(Eigen::MatrixXd theta) {
  // get the mean and variance for each dose group
  Eigen::MatrixXd mu = mean(theta); 
  Eigen::MatrixXd var = variance(theta); 
  
  Eigen::MatrixXd returnV = Y.col(0); // make a log likelihood value
  
  
  // log-likelihood formula in the BMDS documentation - numerically the same as above. 
  if (sufficient_statistics){
    returnV = 		Y.col(2).array()*log(1/sqrt(2.0*M_PI))
    -(Y.col(2).array() / 2.0)*log(var.array()) - (1 / (2.0*var.array()))*((Y.col(2).array() - 1)*pow(Y.col(1).array(), 2) +
                                                                            Y.col(2).array()*pow(Y.col(0).array() - mu.array(), 2));
    return -returnV.sum(); 
  }else{
    Eigen::MatrixXd sqerr = pow(Y.col(0).array()-mu.array(),2); 
    returnV = -(0.5)*log(2*M_PI*var.array())-(1/(2.0*var.array()))*sqerr.array(); 
    return -returnV.sum();

  }
  return 0.0; 
}