#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#ifndef GRADIENT_FUNCTIONS_H
#define GRADIENT_FUNCTIONS_H

template <class PR> 
Eigen::MatrixXd X_logPrior( Eigen::MatrixXd theta, Eigen::MatrixXd p){
      PR prior(p); 
      return prior.log_prior(theta); 
}

void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd,void*)> math_func); 



template <class LL> 
void xgrad2(Eigen::MatrixXd v, double *g, LL *data, Eigen::MatrixXd dose){
    
    Eigen::VectorXd h(v.rows());
    double mpres = pow(1.0e-16, 0.33333);
    double x, temp;
    double derEst;
    Eigen::MatrixXd f1;  Eigen::MatrixXd f2;
    Eigen::MatrixXd  hvector = v;
    
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
    
    /*find the gradients for each of the variables in the likelihood*/
    for (int i = 0; i < v.rows(); i++) {
      /*perform a finite difference calculation on the specific derivative*/
      x = v(i, 0);
      // add h
      hvector(i, 0) = x + h[i];
      f1 = data->mmean(hvector,dose); 
      // subtract h
      hvector(i, 0) = x - h[i];
      f2 = data->mmean(hvector,dose);
      // estimate the derivative
      g[i]  = (f1(0,0) - f2(0,0)) / (2.0*h[i]);
      
      hvector(i, 0) = x;
    }
    return;
}
  
  
template <class LL> 
void xgrad(Eigen::MatrixXd v, double *g, LL *data, Eigen::MatrixXd dose){

  Eigen::VectorXd h(v.rows());
  double mpres = pow(1.0e-16, 0.33333);
  double x, temp;
  double derEst;
  Eigen::MatrixXd f1;  Eigen::MatrixXd f2;
  Eigen::MatrixXd  hvector = v;
  
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
 
  /*find the gradients for each of the variables in the likelihood*/
  for (int i = 0; i < v.rows(); i++) {
    /*perform a finite difference calculation on the specific derivative*/
    x = v(i, 0);
    // add h
    hvector(i, 0) = x + h[i];
    f1 = data->mean(hvector,dose); 
    // subtract h
    hvector(i, 0) = x - h[i];
    f2 = data->mean(hvector,dose);
    // estimate the derivative
    g[i]  = (f1(0,0) - f2(0,0)) / (2.0*h[i]);
   
    hvector(i, 0) = x;
  }
  return;
}


#endif
