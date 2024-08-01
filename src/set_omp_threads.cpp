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
// output: none
// [[Rcpp::export(".set_threads")]]
void set_threads(int num_threads) {
  #ifndef NO_OMP
  if(omp_get_max_threads() > 1){
    if (num_threads > omp_get_num_threads()) {
      Seeder* s = Seeder::getInstance();
      omp_set_num_threads(num_threads);
      s->reset_max_threads(num_threads);
    }
  }else{
    Rcout << "OMP will not be used for parallelization.";
  }
  #else
  Rcpp::Rcout << "OpenMP not supported on this architecure." << std::endl;
  #endif
}
