#include "seeder.h"
#include <nlopt.hpp>
Seeder *Seeder::instance = nullptr;

std::mutex Seeder::instanceMutex;
thread_local gsl_rng *Seeder::r = nullptr;
thread_local int Seeder::currentSeed = -1;
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