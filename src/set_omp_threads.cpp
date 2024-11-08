#ifndef NO_OMP
  #include <omp.h>
#endif

#include <Rcpp.h>
#include "seeder.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
// function: set_threads
// purpose: takes an integer value and set the number of
// openmp threads to that value
// output: integer representing if openMP is supported and set correctely
// [[Rcpp::export(".set_threads")]]
int set_threads(int num_threads) {
  #ifndef NO_OMP
  if(omp_get_max_threads() > 1){
    if (num_threads > omp_get_num_threads()) {
      { omp_set_num_threads(num_threads); }
      // Rcpp::Rcout << "OpenMP threads set to " << num_threads << std::endl;
      return 1;
    }
  }else{
    return 0;
    // Rcpp::Rcout << "OMP will not be used for parallelization.";
  }
  #else
  return -1;
  // Rcpp::Rcout << "OpenMP not supported on this architecure." << std::endl;
  #endif
  return -1;
}
