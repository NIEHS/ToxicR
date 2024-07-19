#include "omp.h"

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
  if (num_threads > omp_get_num_threads()) {
    Seeder* s = Seeder::getInstance();
    omp_set_num_threads(num_threads);
    s->reset_max_threads(num_threads);
  }
}
