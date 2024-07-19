#include "seeder.h"

Seeder* Seeder::instance = nullptr;
std::mutex Seeder::instanceMutex;