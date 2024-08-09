// File: statmodel.h
// Purpose: Create a basic statistical model class
//		   that is used for any general analysis, and can be
//          rather flexible.
// Creator: Matt Wheeler
// Date   : 12/18/2017
// Changes:
//   Matt Wheeler
//   Date 11/06/2018 Added returnX() function to the statModel class
//
//
#include "stdafx.h"
#include <cmath>

#define OPTIM_NO_FLAGS 0
#define OPTIM_USE_GENETIC 1
#define OPTIM_USE_SUBPLX 2
#define OPTIM_USE_BIG_GENETIC 4
#define OPTIM_ALL_FLAGS 7

#ifdef R_COMPILATION
// necessary things to run in R
// necessary things to run in R
// #ifdef ToxicR_DEBUG
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wignored-attributes"
// #include <RcppEigen.h>
// #pragma GCC diagnostic pop
// #else
#include <RcppEigen.h>
// #endif
#include <RcppGSL.h>
#else
#include <Eigen/Dense>

#endif

#include "IDPrior.h"
#include "binomModels.h"
#include "log_likelihoods.h"
#include "seeder.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <nlopt.hpp>

#include <iostream>
#pragma once
#ifndef statmodH
#define statmodH

using namespace std;
struct optimizationResult {
  nlopt::result result;
  double functionV;
  Eigen::MatrixXd max_parms;

  optimizationResult &operator=(const optimizationResult &M) {
    result = M.result;
    functionV = M.functionV;
    max_parms = M.max_parms;
    return *this;
  }
};

template <class LL, class PR> class statModel {
public:
  statModel(LL t_L, PR t_PR, std::vector<bool> b_fixed,
            std::vector<double> d_fixed)
      : log_likelihood(t_L), prior_model(t_PR), isFixed(b_fixed),
        fixedV(d_fixed) {

    if (isFixed.size() != fixedV.size()) {
      throw std::runtime_error(std::string(
          "Statistical Model: Fixed parameter constraints are same size"));
    }

    if (isFixed.size() != (unsigned)log_likelihood.nParms()) {
      throw std::runtime_error(std::string(
          "Statistical Model: Fixed number of parameter constraints not equal "
          "to number of parameters in likelihood model."));
    }
  };

  /////////////////////////////////////////////////////////////////////
  // compute the variance covariance matrix as well as the
  // gradient using finite derivatives.
  Eigen::MatrixXd varMatrix(Eigen::MatrixXd theta);
  Eigen::MatrixXd gradient(Eigen::MatrixXd v);
  //////////////////////////////////////////////////////////////////////

  double negPenLike(Eigen::MatrixXd x, bool bound_check = false) {
    ///////////////////////////////////////////
    for (size_t i = 0; i < isFixed.size(); i++) {
      if (isFixed[i]) {
        x(i, 0) = fixedV[i];
      }
    }
    ///////////////////////////////////////////
    double a = log_likelihood.negLogLikelihood(x);
    double b = prior_model.neg_log_prior(x, bound_check);
    return a + b; // log_likelihood.negLogLikelihood(x) +
                  // prior_model.neg_log_prior(x);
  };

  int nParms() { return log_likelihood.nParms(); };

  Eigen::MatrixXd startValue() { return prior_model.prior_mean(); }

  Eigen::MatrixXd parmLB() { return prior_model.lowerBounds(); }

  Eigen::MatrixXd parmUB() { return prior_model.upperBounds(); }
  // functions to get and set the values of the central estimate of the
  // statistical model
  virtual void setEST(const Eigen::MatrixXd x) {
    Eigen::MatrixXd thetaX = x;
    // worry about the fixed parameters
    for (size_t i = 0; i < isFixed.size(); i++) {
      if (isFixed[i]) {
        thetaX(i, 0) = fixedV[i];
      }
    }
    theta = thetaX;
  };

  virtual Eigen::MatrixXd getEST() {
    Eigen::MatrixXd thetaX = theta;
    // worry about the fixed parameters
    for (size_t i = 0; i < isFixed.size(); i++) {
      if (isFixed[i]) {
        thetaX(i, 0) = fixedV[i];
      }
    }
    return thetaX;
  };

  Eigen::MatrixXd returnX() { return log_likelihood.returnX(); }

public:
  // A stat model has a Log Likelihood and
  // A prior model over the parameters
  // Note: The prior over the parameters does not have to be an actual
  // prior. It can merely place bounds (e.g., box, equality) on functions
  // of the parameters.
  LL log_likelihood;
  PR prior_model;
  // The parameter vector that the log_likelihood and the prior model
  // both depend upon
  std::vector<bool> isFixed;
  std::vector<double> fixedV;
  Eigen::MatrixXd theta;
};

/////////////////////////////////////////////////////////////////
// Function: gradient(Eigen::MatrixXd v)
// Purpose : Performs a finite difference evaluation of the
//           gradient
// Input   :
//          Eigen::MatrixXd v - vector of parameters to model
//
// Output  :
//          double * g      - vector of gradients cooresponding to the
//                            parameters v
///////////////////////////////////////////////////////////////////
template <class LL, class PR>
Eigen::MatrixXd statModel<LL, PR>::gradient(Eigen::MatrixXd v) {

  Eigen::VectorXd h(log_likelihood.nParms());
  double mpres = pow(1.0e-16, 0.5);
  double x, temp;

  double f1 = 0.0;
  double f2 = 0.0;
  Eigen::MatrixXd hvector = v;
  Eigen::MatrixXd g(log_likelihood.nParms(), 1);

  for (int i = 0; i < log_likelihood.nParms(); i++) {
    x = v(i, 0);
    if (fabs(x) > DBL_EPSILON) {
      h[i] = mpres * (fabs(x));
      temp = x + h[i];
      h[i] = temp - x;
    } else {
      h[i] = mpres;
    }
  }
  /*find the gradients for each of the variables in the likelihood*/
  for (int i = 0; i < log_likelihood.nParms(); i++) {
    /*perform a finite difference calculation on the specific derivative*/
    x = v(i, 0);
    // add h
    hvector(i, 0) = x + h[i];
    f1 = negPenLike(hvector);
    // subtract h
    hvector(i, 0) = x - h[i];
    f2 = negPenLike(hvector);
    // estimate the derivative
    g(i, 0) = (f1 - f2) / (2.0 * h[i]);
    hvector(i, 0) = x;
  }
  return g;
}

////////////////////////////////////////////////////////////////////
// Function: varMatrix(Eigen::MatrixXd theta)
// Purpose:  Computes the inverse hessian of the penalized likelihood
//			 given the parameter theta.
// Input  :  Eigen::MatrixXd theta kx1 parameter value equal to the
//           number of parameters in the "statistical model"
// Output :  Eigen::MatrixXd k x k symmetric positive definite matrix
//			 Inv Hessian of the penalized likelihood.
/////////////////////////////////////////////////////////////////////
template <class LL, class PR>
Eigen::MatrixXd statModel<LL, PR>::varMatrix(Eigen::MatrixXd theta) {

  Eigen::MatrixXd m(log_likelihood.nParms(), log_likelihood.nParms());
  m.setZero();
  double temp;
  double h = 1e-10;
  double mpres = pow(1.0e-16, 0.333333);

  Eigen::MatrixXd tempVect(log_likelihood.nParms(), 1);
  tempVect.setZero();
  double hi, hj;

  for (int i = 0; i < log_likelihood.nParms(); i++) {
    for (int j = 0; j < log_likelihood.nParms(); j++) {

      double x = theta(i, 0); // compute the hessian for parameter i-j
      if (fabs(x) > DBL_EPSILON) {
        h = mpres * (fabs(x));
        temp = x + h;
        hi = temp - x;
      } else {
        hi = mpres;
      }

      x = theta(j, 0);
      if (fabs(x) > DBL_EPSILON) {
        h = mpres * (fabs(x));
        temp = x + h;
        hj = temp - x;
      } else {
        hj = mpres;
      }

      temp = 0.0;

      if (i == j) {
        // on diagonal variance estimate
        // eg the 2nd partial derivative
        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) + 2 * hi;
        temp += -1.0 * negPenLike(tempVect, true);
        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) + hi;
        temp += 16.0 * negPenLike(tempVect, true);
        temp += -30.0 * negPenLike(theta, true);
        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) - hi;
        temp += 16.0 * negPenLike(tempVect, true);

        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) - 2 * hi;
        temp += -1.0 * negPenLike(tempVect, true);
        temp = temp / (12 * hi * hi);
        m(i, i) = temp;
        // if nan, shift up off boundary
        if (std::isnan(temp)) {
          temp = 0.0;
          tempVect = theta;
          tempVect(i, 0) += 2.0 * hi;
          tempVect(i, 0) = tempVect(i, 0) + 2 * hi;
          temp += -1.0 * negPenLike(tempVect, true);
          tempVect = theta;
          tempVect(i, 0) += 2.0 * hi;
          tempVect(i, 0) = tempVect(i, 0) + hi;
          temp += 16.0 * negPenLike(tempVect, true);
          temp += -30.0 * negPenLike(theta, true);
          tempVect = theta;
          tempVect(i, 0) += 2.0 * hi;
          tempVect(i, 0) = tempVect(i, 0) - hi;
          temp += 16.0 * negPenLike(tempVect, true);
          tempVect = theta;
          tempVect(i, 0) += 2.0 * hi;
          tempVect(i, 0) = tempVect(i, 0) - 2 * hi;
          temp += -1.0 * negPenLike(tempVect, true);
          temp = temp / (12 * hi * hi);
          m(i, i) = temp;
        }
      } else {

        // off diagonal variance estimate
        // eg the first partial of x wrt the
        // first partial  of y

        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) + hi;
        tempVect(j, 0) = tempVect(j, 0) + hj;
        temp += negPenLike(tempVect, true);

        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) + hi;
        tempVect(j, 0) = tempVect(j, 0) - hj;
        temp += -1.0 * negPenLike(tempVect, true);

        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) - hi;
        tempVect(j, 0) = tempVect(j, 0) + hj;
        temp += -1.0 * negPenLike(tempVect, true);

        tempVect = theta;
        tempVect(i, 0) = tempVect(i, 0) - hi;
        tempVect(j, 0) = tempVect(j, 0) - hj;
        temp += negPenLike(tempVect, true);
        temp = temp / (4 * hi * hj);
        m(i, j) = temp;
        // if nan, shift up off boundary
        if (std::isnan(temp)) {
          temp = 0.0;
          tempVect = theta;
          tempVect(i, 0) += hi;
          tempVect(j, 0) += hj;
          tempVect(i, 0) = tempVect(i, 0) + hi;
          tempVect(j, 0) = tempVect(j, 0) + hj;
          temp += negPenLike(tempVect, true);

          tempVect = theta;
          tempVect(i, 0) += hi;
          tempVect(j, 0) += hj;
          tempVect(i, 0) = tempVect(i, 0) + hi;
          tempVect(j, 0) = tempVect(j, 0) - hj;
          temp += -1.0 * negPenLike(tempVect, true);

          tempVect = theta;
          tempVect(i, 0) += hi;
          tempVect(j, 0) += hj;
          tempVect(i, 0) = tempVect(i, 0) - hi;
          tempVect(j, 0) = tempVect(j, 0) + hj;
          temp += -1.0 * negPenLike(tempVect, true);

          tempVect = theta;
          tempVect(i, 0) += hi;
          tempVect(j, 0) += hj;
          tempVect(i, 0) = tempVect(i, 0) - hi;
          tempVect(j, 0) = tempVect(j, 0) - hj;
          temp += negPenLike(tempVect, true);
          temp = temp / (4 * hi * hj);
          m(i, j) = temp;
        }
      }
    }
  }
  // replace NaN by 0
  m = m.unaryExpr([](double xh) { return std::isnan(xh) ? 0.0 : xh; });

  // m = m.inverse();
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(m);
  auto rank = lu_decomp.rank();
  // matrix is less than full rank
  // so we add an epsolon error to fix it
  if (rank < m.rows()) {
    m = m + 1e-4 * Eigen::MatrixXd::Identity(m.rows(), m.cols());
  }
  // return m.inverse();
  m = m.inverse();
  // fix the covariance element if it was zero (and 1e-4 was added)
  m = m.unaryExpr([](double x) { return (x == 1e4) ? 0.0 : x; });

  return m;
}

/////////////////////////////////////////////////////////////////
// Function: neg_pen_likelihood(unsigned n,
//								const double *b,
//								double *grad,
//								void *data)
// Purpose: Used to estimate the negative penalized likelihood from the statMod
// class
//          in the proper NLOPT format
// Input:
//          unsigned n: The length of the vector of parameters to the optimizer
//          const double *b   : A pointer to the place in memory the vector
//          resides double *gradient  : The gradient of the function at b, NULL
//          if unused.
//                              It is currently unused.
//          void    *data     : Extra data needed. In this case, it is a
//          statModel<LL,PR> object,
//							   which is used to
//compute the negative penalized likelihood
//////////////////////////////////////////////////////////////////
template <class LL, class PR>
double neg_pen_likelihood(unsigned n, const double *b, double *grad,
                          void *data) {
  statModel<LL, PR> *model = (statModel<LL, PR> *)data;
  Eigen::MatrixXd theta(n, 1);
  for (unsigned int i = 0; i < n; i++)
    theta(i, 0) = b[i];

  // if an optimizer is called that needs a gradient we calculate
  // it using the finite difference gradient
  if (grad) {
    Eigen::MatrixXd mgrad = model->gradient(theta);
    for (int i = 0; i < model->nParms(); i++) {
      grad[i] = mgrad(i, 0);
    }
  }

  return model->negPenLike(theta);
}

////////////////////////////////////////////////////////////
// Function: startValue_F(statModel<LL, PR>  *M)
// Purpose : Given a statistical model M it produces a starting
//           value based upon one iteration of a Fisher scoring algorithm.
//           As there is no reason to assume the single iteration will respect
//			the parameter bounds, the routine checks if the
// parameter bounds 		    are violated; if any are, it moves in the
// direction of the change
//           in half power increments until the bounds are not violated. This
//           value
//	        is returned.
// Input    : statModel<LL,PR> *M - Pointer to a statistical value whose central
//			 tendency values in PR respect the constraints of the
// model. Output   : A starting value to the optimizer that *should* allow the
// optimizer to
//          : proceed quickly to the optimum.
/////////////////////////////////////////////////////////////

template <class LL, class PR>
std::vector<double> startValue_F(statModel<LL, PR> *M, Eigen::MatrixXd startV,
                                 std::vector<double> lb, std::vector<double> ub,
                                 bool isBig = true) {

  std::vector<double> x(M->nParms());

  int NI;
  if (isBig) {
    NI = 1000; // size of the initial population
  } else {
    NI = 500;
  }
  std::vector<double> llist(
      NI + 1, std::numeric_limits<double>::infinity()); // List of the
                                                        // likelihood values;
  /*for(int j=0; j< llist.size(); j++){
    llist[j] = std::numeric_limits<double>::infinity(); // initialize everything
  to infinity
  }*/
  std::vector<Eigen::MatrixXd> population(
      NI + 1); // list of the population parameters

  double test_l;
  // make sure start value is within our bounds
  for (unsigned int i = 0; i < lb.size(); i++) {
    if (startV(i, 0) < lb[i] || startV(i, 0) > ub[i]) {
      startV(i, 0) = lb[i];
    }
  }
  Eigen::MatrixXd test = startV;
  // test = M->startValue();
  Seeder *seeder = Seeder::getInstance();

  population[NI] = startV;
  llist[NI] = M->negPenLike(test);

  //
  // create the initial population of size (NI) random starting points for the
  // genetic algorithm double initial_temp;
  for (int i = 0; i < NI; i++) {
    // generate new value to be within the specified bounds
    for (int j = 0; j < M->nParms(); j++) {
      test(j, 0) =
          startV(j, 0) + seeder->get_ran_flat(); // random number in the bounds

      if (test(j, 0) > ub[j]) {
        test(j, 0) = ub[j];
      }
      if (test(j, 0) < lb[j]) {
        test(j, 0) = lb[j];
      }
    }
    test_l = M->negPenLike(test);
    bool break_loop = false;
    // put the new value in sorted order based upon likelihood
    // score
    for (int j = 0; !break_loop && j < NI; j++) {
      if (test_l < llist[j]) { // this is the first occurance
        std::vector<double>::iterator it_l = llist.begin();
        std::vector<Eigen::MatrixXd>::iterator it_pop = population.begin();

        std::advance(it_l, j);
        std::advance(it_pop, j);
        llist.insert(it_l, test_l);
        population.insert(it_pop, test);
        break_loop = true;
      }
    }
  }
  // look for bad population entries
  for (int i = population.size() - 1; i > 1; --i) {
    if (population[i].size() == 0) {
      population.erase(population.begin() + i);
      i = population.size(); // removed the value
                             // start over
    }
  }

  if (population.size() <= 25) {
    // couln't find a good starting point return the starting value
    // and pray
    for (int i = 0; i < M->nParms(); i++)
      x[i] = startV(i, 0);
    return x;
  }
  //
  // Now do the Genetic algoritm thing.
  // first trim the population to allow only
  // the fittest to procrate
  std::vector<double>::iterator it_l = llist.begin();
  std::vector<Eigen::MatrixXd>::iterator it_pop = population.begin();
  int tmp = std::min<int>(population.size(), 175);
  std::advance(it_l, tmp);
  std::advance(it_pop, tmp);
  llist.erase(it_l, llist.end());
  population.erase(it_pop, population.end());

  int ngenerations;
  ;
  int ntourny;
  int tourny_size;
  if (isBig) {
    ngenerations = 600;
    ntourny = 30;
    tourny_size = 40;
  } else {
    ngenerations = 450;
    ntourny = 20;
    tourny_size = 20;
  }

  for (int xx = 0; xx < ngenerations; xx++) {

    std::vector<double> tourny_winners;
    std::vector<Eigen::MatrixXd> tourny_vals;

    for (int ay = 0; ay < ntourny; ay++) {

      std::vector<double> cur_tourny_nll(
          tourny_size, std::numeric_limits<double>::infinity());
      std::vector<Eigen::MatrixXd> cur_tourny_parms(tourny_size, population[0]);
      // make sure there are no entries in the current tourny.

      // select the elements from the population
      for (int z = 0; z < tourny_size; z++) {
        // in each generation there is a tournament
        // find the best individual out of tourny_size this individual is:
        // randomly mutatied and differentially evolved based upon the given
        // individuals in the tourny.
        int sel = (int)(population.size() *
                        seeder->get_uniform()); // choose which element in the
                                                // population
        cur_tourny_nll[z] = it_l[sel];
        cur_tourny_parms[z] = population[sel];
      }

      // find the best
      double best_nll = cur_tourny_nll[0];
      Eigen::MatrixXd best_parm = cur_tourny_parms[0];

      for (int z = 1; z < tourny_size; z++) {
        if (cur_tourny_nll[z] < best_nll) {
          best_nll = cur_tourny_nll[z];
          best_parm = cur_tourny_parms[0];
        }
      }
      // the best is the zero element
      // randomly select another element to find the diference
      bool correctBounds = true;

      int idx = (int)(cur_tourny_parms.size() - 1) * seeder->get_uniform() + 1;
      Eigen::MatrixXd temp_delta = best_parm - cur_tourny_parms[idx];
      // Create a new child as a mix between the best and some other
      // value.
      Eigen::MatrixXd child =
          best_parm + 0.8 * temp_delta * (2 * seeder->get_uniform() - 1);
      correctBounds = true;

      for (int iii = 0; iii < M->nParms(); iii++) {
        // perterb the individual values in the child
        child(iii, 0) = child(iii, 0) + 0.2 * abs(child(iii, 0)) *
                                            (2 * seeder->get_uniform() - 1);
        if (lb[iii] > child(iii, 0) || ub[iii] < child(iii, 0)) {
          correctBounds = false;
          break;
        }
      }

      if (correctBounds) {
        test_l = M->negPenLike(child);
      } else {
        test_l = std::numeric_limits<double>::infinity();
      }

      // put this new child into the population
      bool break_loop = false;
      int S = population.size();

      for (int j = 0; !break_loop && j < S; j++) {
        if (test_l < llist[j]) { // this is the first occurance
          std::vector<double>::iterator it_l = llist.begin();
          std::vector<Eigen::MatrixXd>::iterator it_pop = population.begin();

          std::advance(it_l, j);
          std::advance(it_pop, j);
          llist.insert(it_l, test_l);
          population.insert(it_pop, child);
          break_loop = true;
        }
      }
    }

    it_l = llist.begin();
    it_pop = population.begin();

    if (llist.size() >= 100)
      std::advance(it_l, 100);
    else
      std::advance(it_l, llist.size());

    if (population.size() >= 100)
      std::advance(it_pop, 100);
    else
      std::advance(it_pop, population.size());

    llist.erase(it_l, llist.end());
    population.erase(it_pop, population.end());
  }

  if (population.size() > 0) {
    test = population[0]; // the fittest is our starting value
  }

  double t1 = M->negPenLike(test);
  double t2 = M->negPenLike(startV);

  if (t2 < t1) { // the random search was no better than the first value.
    test = startV;
  }
  bool found_nans = false;
  for (int i = 0; i < M->nParms(); i++) {
    if (isnan(test(i, 0))) {
      found_nans = true;
    }
  }
  if (found_nans) { // something really messed up revert to initial
                    // starting values
    test = startV;
  }

  for (int i = 0; i < M->nParms(); i++)
    x[i] = test(i, 0);

  for (int i = 0; i < M->nParms(); i++) {
    if (!isnormal(x[i])) {
      x[i] = 0;
    }
  }
  return x;
}

////////////////////////////////////////////////////////////////////
// Function: findMAP(statMod<LL,PR> *M)
// Find the maximum a posteriori estimate given a statistical
// model
// Input : statMod<LL,PR> *M - A given statistical model;
//  note by default all optimization flags are on.
// Output: statMod<LL,PR> *M - The model with it's MAP parameter set.
template <class LL, class PR>
optimizationResult findMAP(statModel<LL, PR> *M, Eigen::MatrixXd startV,
                           unsigned int flags = OPTIM_USE_GENETIC |
                                                OPTIM_USE_SUBPLX) {
  optimizationResult oR;
  Eigen::MatrixXd temp_data = M->parmLB();
  std::vector<double> lb(M->nParms());

  for (int i = 0; i < M->nParms(); i++)
    lb[i] = temp_data(i, 0);
  temp_data = M->parmUB();

  std::vector<double> ub(M->nParms());

  for (int i = 0; i < M->nParms(); i++)
    ub[i] = temp_data(i, 0);
  std::vector<double> x(startV.rows());
  if (OPTIM_USE_GENETIC & flags) {
    bool op_size = (OPTIM_USE_BIG_GENETIC & flags);
    try {

      x = startValue_F(M, startV, lb, ub, op_size);

    } catch (...) {
    }
  } else {
    for (unsigned int i = 0; i < x.size(); i++) {
      x[i] = startV(i, 0);
    }
  }

  // int yy = x.size();

  for (int i = 0; i < M->nParms(); i++) {
    if (!isnormal(x[i])) {
      x[i] = 0;
    }
  }

  double minf;
  nlopt::result result = nlopt::FAILURE;

  std::vector<double> init(x.size());
  // for (int i = 0; i < x.size(); i++) init[i] = 1e-4;
  //  set up a bunch of differnt plausible optimizers in case of failure
  //  the first one is mainly to get a better idea of a starting value
  //  though it often converges to the optimum.

  nlopt::opt opt1(nlopt::LN_SBPLX, M->nParms());
  nlopt::opt opt3(nlopt::LD_LBFGS, M->nParms());
  nlopt::opt opt2(nlopt::LN_BOBYQA, M->nParms());

  nlopt::opt opt4(nlopt::LN_COBYLA, M->nParms());
  nlopt::opt opt5(nlopt::LD_SLSQP, M->nParms());

  nlopt::opt *opt_ptr;

  int opt_iter;
  // look at 5 optimization algorithms :-)
  int start_iter = 0; //(OPTIM_USE_SUBPLX & flags)?0:1;

  for (opt_iter = start_iter; opt_iter <= 4; opt_iter++) {

    // Ensure that starting values are within bounds
    for (int i = 0; i < M->nParms(); i++) {
      double temp = x[i];
      if (temp < lb[i])
        temp = lb[i];
      else if (temp > ub[i])
        temp = ub[i];
      x[i] = temp;
    } // end for

    switch (opt_iter) {
    case 0:
      opt_ptr = &opt1;
      opt_ptr->set_maxeval(1200);
      break;
    case 1:
      opt_ptr = &opt2;
      opt_ptr->set_maxeval(5000);
      break;
    case 2:
      opt_ptr = &opt3;
      opt_ptr->set_maxeval(5000);
      break;
    case 3:
      opt_ptr = &opt4;
      opt_ptr->set_maxeval(5000);
    default:
      opt_ptr = &opt5;
      opt_ptr->set_maxeval(5000);
      break;
    }

    opt_ptr->set_lower_bounds(lb);
    opt_ptr->set_upper_bounds(ub);
    // opt_ptr->set_ftol_rel(1e-8);
    opt_ptr->set_xtol_rel(1e-9);
    // opt_ptr->set_initial_step(1e-3);
    opt_ptr->set_min_objective(neg_pen_likelihood<LL, PR>, M);

    ////////////////////////////////////////////////
    //////////////////////////////////////////////
    // set the start distance to a size
    // nlopt's default options suck
    ///////////////////////////////////////////////

    try {
      result = opt_ptr->optimize(x, minf);
      // note even if it doesn't converge the x will be updated
      // to the best value the optimizer ever had, which will allow
      // the next optimizer to carry on.
      // Exit the loop if good result and not first try.
      // lco: should change to a break statement since "10"
      // could be legitimate iteration in the future

      if (opt_iter >= 1 && result > 0 && result < 5) {

        opt_iter = 10; // if it made it here it will break the loop
      }

    } // try
    catch (const std::invalid_argument &exc) {
      DEBUG_LOG(file, "opt_iter= " << opt_iter
                                   << ", error: invalid arg: " << exc.what());

    } // catch
    catch (nlopt::roundoff_limited &exec) {
      DEBUG_LOG(file, "opt_iter= " << opt_iter << ", error: roundoff_limited");
      //  cout << "bogo" << endl;
    } // catch
    catch (nlopt::forced_stop &exec) {
      DEBUG_LOG(file, "opt_iter= " << opt_iter << ", error: forced_stop");
      //  cout << "there" << endl;
    } // catch
    catch (const std::exception &exc) {
      DEBUG_LOG(file,
                "opt_iter= " << opt_iter << ", general error: " << exc.what());
      // cout << "???" << endl;
    } catch (...) {

    } // catch

    DEBUG_CLOSE_LOG(file);
  } // for opt_iter

  Eigen::Map<Eigen::MatrixXd> d(x.data(), M->nParms(), 1);
  oR.result = result;
  oR.functionV = minf;
  oR.max_parms = d;
  // cout << "Result Opt:" << result << "\n---\n"<< endl << d ;

  M->setEST(d);

  if (result < 0) {
    // cerr << __FUNCTION__ << " at line: " << __LINE__ << " result= " << result
    // << endl;
  }

  return oR;
}

////////////////////////////////////////////////////////////////////
// Function: findMAP(statMod<LL,PR> *M)
// Find the maximum a posteriori estimate given a statistical
// model
// Input : statMod<LL,PR> *M - A given statistical model;
// Output: statMod<LL,PR> *M - The model with it's MAP parameter set.
template <class LL, class PR> optimizationResult findMAP(statModel<LL, PR> *M) {

  return findMAP<LL, PR>(M, M->startValue());
}

template <class LL, class PR>
optimizationResult findMAP(statModel<LL, PR> *M, signed int flags) {
  return findMAP<LL, PR>(M, M->startValue(), flags);
}

#endif
