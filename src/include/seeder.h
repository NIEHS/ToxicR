#define STRICT_R_HEADERS
#pragma once
#ifndef SEEDER
#define SEEDER

#ifndef NO_OMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <mutex>
#include <thread>
#include <vector>

class Seeder {
private:
  gsl_rng *r = nullptr;
  static Seeder *instance;
  static std::mutex instanceMutex;
  int max_threads;
  std::mutex seedMutex;
  Seeder() {}

public:
  int currentSeed;
  static thread_local gsl_rng *r_local;
  std::vector<gsl_rng *> rngs;
  const gsl_rng_type *T;
  Seeder(Seeder const &) = delete;
  Seeder &operator=(Seeder const &) = delete;

  ~Seeder() {
    for (int i = 0; i < max_threads; i++) {
      gsl_rng_free(rngs[i]);
    }
  }
#ifndef NO_OMP // OpenMP - Multi-threaded seeder
  static Seeder *getInstance() {
    std::lock_guard<std::mutex> lock(instanceMutex);
    if (!instance) {
      instance = new Seeder();
      omp_set_dynamic(0); // Disable dynamic threads
      omp_set_nested(0);  // Disable nested parallelism
      instance->max_threads = omp_get_num_threads();
      instance->T = gsl_rng_mt19937;
      instance->currentSeed = 0;
      instance->rngs.resize(instance->max_threads, nullptr);
    }
    return instance;
  }

  void setSeed(int seed) {
    if (seed < 0) {
      Rcpp::stop("Error: Seed must be a positive integer.");
    }
    currentSeed = seed;
    if (r_local) {
      gsl_rng_free(r_local);
    }

    r_local = gsl_rng_alloc(T);
    gsl_rng_set(r_local, seed);
  }

  double get_uniform() {
    if (!r_local) {
      r_local = gsl_rng_alloc(T);
      gsl_rng_set(r_local, currentSeed);
    }
    return gsl_rng_uniform(r_local);
  }

  double get_gaussian_ziggurat() {

    if (!r_local) {
      r_local = gsl_rng_alloc(T);
      gsl_rng_set(r_local, currentSeed);
    }
    return gsl_ran_gaussian_ziggurat(r_local, 1.0);
  }

  double get_ran_flat() {

    if (!r_local) {
      r_local = gsl_rng_alloc(T);
      gsl_rng_set(r_local, currentSeed);
    }
    return gsl_ran_flat(r_local, -1, 1);
  }
#else // Non-OpenMP - single-threaded
  static Seeder *getInstance() {
    if (!instance) {
      instance = new Seeder();
      instance->max_threads = 1;
      instance->T = gsl_rng_mt19937;
      instance->currentSeed = 0;
      instance->rngs.resize(instance->max_threads, nullptr);
      gsl_rng *r_local = gsl_rng_alloc(T);
      gsl_rng_set(r_local, currentSeed);
      rngs[thread_num] = r_local;
    }
    return instance;
  }

  void setSeed(int seed) {
    if (seed < 0) {
      Rcpp::stop("Error: Seed must be a positive integer.");
    }
    currentSeed = seed;
    if (r_local) {
      gsl_rng_free(r_local);
    }
    r_local = gsl_rng_alloc(T);
    gsl_rng_set(r_local, seed);
  }

  double get_uniform() { return gsl_rng_uniform(r_local); }

  double get_gaussian_ziggurat() {
    return gsl_ran_gaussian_ziggurat(r_local, 1.0);
  }

  double get_ran_flat() { return gsl_ran_flat(r_local, -1, 1); }

#endif // NO_OMP
};
#endif // SEEDER