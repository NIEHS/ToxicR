#include "seeder.h"
#include <nlopt.hpp>

// [[Rcpp::depends(RcppGSL)]]
// function: setseedGSL
// purpose: takes an integer value and sets the GSL seed
// output: none
// [[Rcpp::export(".setseedGSL")]]
void setseedGSL(int s) {
  Seeder *seeder = Seeder::getInstance();
  seeder->setSeed(s);
  nlopt_srand(s);
  return;
}