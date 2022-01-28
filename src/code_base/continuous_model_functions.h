
#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

Eigen::MatrixXd standard_data(Eigen::MatrixXd Y,Eigen::MatrixXd X,double *zero_mean);
