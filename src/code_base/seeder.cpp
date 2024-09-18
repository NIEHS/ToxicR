#include "seeder.h"

thread_local gsl_rng* Seeder::r_local = nullptr;
Seeder* Seeder::instance = nullptr;
std::mutex Seeder::instanceMutex;