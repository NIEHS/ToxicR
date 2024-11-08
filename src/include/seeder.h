#define STRICT_R_HEADERS
#pragma once
#ifndef SEEDER
#define SEEDER

// Check for MinGW or MSVC (both support __declspec(thread))
#if defined(_MSC_VER) || defined(__MINGW32__)
#define THREAD_LOCAL __declspec(thread)

// Check for GCC or Clang on non-Windows systems
#elif defined(__GNUC__) || defined(__clang__)
#define THREAD_LOCAL thread_local

#else
#define THREAD_LOCAL // Fallback for unsupported compilers
#endif
#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <mutex>
#include <thread>
#include <vector>
#include <nlopt.hpp>

class Seeder {
private:
  THREAD_LOCAL static gsl_rng *rng;
  const gsl_rng_type *T = gsl_rng_mt19937;
  int currentSeed = 0;

  Seeder() {}
  Seeder(Seeder const &) = delete;
  Seeder &operator=(Seeder const &) = delete;

public:
  static Seeder *getInstance() {
    static Seeder instance;
    return &instance;
  }

  ~Seeder() {
    if (rng) {
      gsl_rng_free(rng);
    }
  }

  void setSeed(int seed) {
    if (seed < 0) {
      Rcpp::stop("Error: Seed must be a positive integer.");
    }
    if (!rng) {
      rng = gsl_rng_alloc(T);
    }
    gsl_rng_set(rng, seed);
    // Rcpp::Rcout << "GSL seed set to: " << seed << std::endl;
    nlopt_srand(seed);
    currentSeed = seed;
  }

  double get_uniform() {
    if (!rng) {
      Rcpp::warning("Error: RNG not initialized.");
      setSeed(currentSeed);
    }
    return gsl_rng_uniform(rng);
  }

  double get_gaussian_ziggurat() {
    if (!rng) {
      Rcpp::warning("Error: RNG not initialized.");
      setSeed(currentSeed);
    }
    return gsl_ran_gaussian_ziggurat(rng, 1.0);
  }

  double get_ran_flat() {
    if (!rng) {
      Rcpp::warning("Error: RNG not initialized.");
      setSeed(currentSeed);
    }
    return gsl_ran_flat(rng, -1, 1);
  }
};
#endif