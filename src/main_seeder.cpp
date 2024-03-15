#include "seeder.h"

Seeder* Seeder::instance = nullptr;

// [[Rcpp::depends(RcppGSL)]]
// function: setseedGSL
// purpose: takes an integer value and sets the GSL seed
// output: none
// [[Rcpp::export(".setseedGSL")]]
void setseedGSL(int s) {
  Seeder* seeder = Seeder::getInstance();
  seeder->setSeed(s);
  return;
}