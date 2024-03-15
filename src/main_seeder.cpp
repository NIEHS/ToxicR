#include "seeder.h"

Seeder* Seeder::instance = nullptr;

// [[Rcpp::export]]
void setseedGSL(int s) {
  Seeder* seeder = Seeder::getInstance();
  seeder->setSeed(s);
  return;
}