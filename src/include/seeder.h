#define STRICT_R_HEADERS
#pragma once
#ifndef SEEDER
#define SEEDER

#include "omp.h"
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
  std::vector<gsl_rng*> rngs;
  const gsl_rng_type * T;
  Seeder(Seeder const &) = delete;
  Seeder &operator=(Seeder const &) = delete;
  static Seeder *getInstance() {
    std::lock_guard<std::mutex> lock(instanceMutex);
    if (!instance) {
      instance = new Seeder();
      instance->max_threads = omp_get_num_threads();
      instance->T = gsl_rng_mt19937;
      instance->currentSeed = 0;
      instance->rngs.reserve(instance->max_threads);
    }

    return instance;
  }

  void reset_max_threads(int threads) {
    if(max_threads < threads) {
      int num_prev_threads = max_threads;
      max_threads = threads;
      
      rngs.reserve(threads);
      #pragma omp parallel for
      for (int i = 0; i < threads; i++) {
        int thread_num = omp_get_thread_num();
        if (thread_num < num_prev_threads) {
          continue;
        }
        thread_local gsl_rng* r_local = gsl_rng_alloc(T);
        gsl_rng_set(r_local, currentSeed);
        rngs[thread_num] = r_local;
      }
    }
  }

  ~Seeder() {
    for(int i = 0; i < max_threads; i++) {
      gsl_rng_free(rngs[i]);
    }
  }

  void setSeed(int seed) {
    if (seed < 0) {
      Rcpp::stop("Error: Seed must be a positive integer.");
    }
    if(rngs.empty()) {
      #pragma omp parallel for
      for (int i = 0; i < max_threads; i++) {
        int thread_num = omp_get_thread_num();
        thread_local gsl_rng* r_local = gsl_rng_alloc(T);
        gsl_rng_set(r_local, seed);
        rngs[thread_num] = r_local;
      }
    }
    else {
      #pragma omp parallel for
      for (int i = 0; i < max_threads; i++) {
        int thread_num = omp_get_thread_num();
        gsl_rng_free(rngs[thread_num]);
        thread_local gsl_rng* r_local = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r_local, seed);
        rngs[thread_num] = r_local;
      }
    }
  }
  
  double get_uniform() {
    double r_val;
    int thread_num = omp_get_thread_num();
    r_val = gsl_rng_uniform(rngs[thread_num]);
    return r_val;
  }

  double get_gaussian_ziggurat() {
    double r_val;
    int thread_num = omp_get_thread_num();
    r_val = gsl_ran_gaussian_ziggurat(rngs[thread_num], 1.0);
    return r_val;
  }

  double get_ran_flat() {
    double r_val;
    int thread_num = omp_get_thread_num();
    r_val = gsl_ran_flat(rngs[thread_num], -1, 1);
    return r_val;
  }
};
#endif