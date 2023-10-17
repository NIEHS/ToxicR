#include <omp.h>

// function: set_threads
// purpose: takes an integer value and set the number of
// openmp threads to that value
// output: none
// [[Rcpp::export(".set_threads")]]
void set_threads(int num_threads) {
    if (num_threads != omp_get_num_threads()) {
        omp_set_num_threads(num_threads);
    }
}