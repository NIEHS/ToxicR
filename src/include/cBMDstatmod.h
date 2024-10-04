// File: cBMDstatmod.h
// Purpose: Creates a basic statistical model class
//          that is a BMD statistical model. This class provides
//          all of the information necessary to compute a benchmark dose
//          for continuous data
// Creator: Matt Wheeler
// Date   : 4/19/2018
// Changes:
//

#include "cmodeldefs.h"
#include "statmod.h"

#pragma once
#ifndef cBMDstatmodH
#define cBMDstatmodH

#include <cmath>
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

#include <list>

#include <iostream>
#include <limits>
#include <nlopt.hpp>

using namespace std;

template <class LL, class PR> // log-likehood, prior model
class cBMDModel : public statModel<LL, PR> {
public:
  cBMDModel(LL t_L, PR t_PR, std::vector<bool> b_fixed,
            std::vector<double> d_fixed, bool I)
      : statModel<LL, PR>(t_L, t_PR, b_fixed, d_fixed), isInc(I){};

  bool isIncreasing() { return isInc; }

  int type_of_profile(contbmd TYPE) {
    return this->log_likelihood.type_of_profile(TYPE);
  };

  int parameter_to_remove(contbmd TYPE) {
    return this->log_likelihood.parameter_to_remove(TYPE);
  };

  cont_model modelling_type() { return this->log_likelihood.mean_type(); }

  ////////////////////////////////////////////////////////////////////
  // BASIC BMD MANIPULATIONS
  double returnBMD(contbmd BMDType, double BMRF, double advP) {
    return returnBMD(this->theta, BMDType, BMRF, advP);
  };
  double returnBMD(Eigen::MatrixXd theta, contbmd BMDType, double BMRF,
                   double advP);
  ////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  double equality_constraint(Eigen::MatrixXd theta, double *grad,
                             contbmd BMDType, double BMD, double BMRF,
                             double advP);

  ////////////////////////////////////////////////////////////////////
  //
  //
  //
  ////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd starting_value(Eigen::MatrixXd theta, contbmd BMDType,
                                 double BMD, double BMRF, bool isInc,
                                 double tail_prob, std::vector<double> lb,
                                 std::vector<double> ub) {

    return this->log_likelihood.starting_value(theta, BMDType, BMD, BMRF, isInc,
                                               tail_prob, lb, ub);
  };

  std::vector<double> bound_fix(std::vector<double> x, contbmd BMDType,
                                double BMRF, double tail_prob, double BMD,
                                bool isInc);

private:
  bool isInc;
};
//////////////////////////
template <class A, class B> struct c_optimInfo {
  cBMDModel<A, B> *sm;
  double cBMD;
  double BMR;
  contbmd BMDType;
  double add_info;
  bool isExtra;
  bool isInc;
};
////////////////////

////////////////////
template <class LL, class PR>
std::vector<double> cBMDModel<LL, PR>::bound_fix(std::vector<double> x,
                                                 contbmd BMDType, double BMRF,
                                                 double tail_prob, double BMD,
                                                 bool isInc) {
  switch (BMDType) {
  case CONTINUOUS_BMD_ABSOLUTE:
    return this->log_likelihood.bmd_start_absolute_clean(x, BMRF, BMD, isInc);
    break;
  case CONTINUOUS_BMD_STD_DEV:

    return this->log_likelihood.bmd_start_stddev_clean(x, BMRF, BMD, isInc);
    break;
  case CONTINUOUS_BMD_REL_DEV:
    return this->log_likelihood.bmd_start_reldev_clean(x, BMRF, BMD, isInc);
    break;
  case CONTINUOUS_BMD_POINT:
    return this->log_likelihood.bmd_start_point_clean(x, BMRF, BMD, isInc);
    break;
  case CONTINUOUS_BMD_EXTRA:
    return this->log_likelihood.bmd_start_extra_clean(x, BMRF, BMD, isInc);
    break;
  case CONTINUOUS_BMD_HYBRID_EXTRA:
    return this->log_likelihood.bmd_start_hybrid_extra_clean(x, BMRF, BMD,
                                                             isInc, tail_prob);
    break;
  case CONTINUOUS_BMD_HYBRID_ADDED:

  default:
    return x;
    break;
  }
}
/*
 *Function: double equality_constraint(Eigen::MatrixXd theta, double *grad,
 *double BMD, double BMRF, double advP); Purpose : Give the current value of the
 *inequality constraint Output: Double that describes the current value of the
 *inequality
 */
template <class LL, class PR>
double cBMDModel<LL, PR>::equality_constraint(Eigen::MatrixXd theta,
                                              double *grad, contbmd BMDType,
                                              double BMD, double BMRF,
                                              double advP) {

  // set all the fixed parameters to their value
  for (size_t i = 0; i < this->isFixed.size(); i++) {
    if (this->isFixed[i]) {
      theta(i, 0) = this->fixedV[i];
    }
  }

  if (grad != NULL) {
    Eigen::MatrixXd GRAD = this->log_likelihood.eqConst_gradient(
        theta, BMDType, BMD, BMRF, true, advP);
    for (int i = 0; i < theta.rows(); i++)
      grad[i] = GRAD(i, 0);
  }
  double returnV = 0.0;
  switch (BMDType) {
  case CONTINUOUS_BMD_ABSOLUTE:
    returnV = this->log_likelihood.bmd_absolute_bound(theta, BMD, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_STD_DEV:
    returnV = this->log_likelihood.bmd_stdev_bound(theta, BMD, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_REL_DEV:
    returnV = this->log_likelihood.bmd_reldev_bound(theta, BMD, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_POINT:
    returnV = this->log_likelihood.bmd_point_bound(theta, BMD, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_EXTRA:
    returnV = this->log_likelihood.bmd_extra_bound(theta, BMD, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_HYBRID_EXTRA:
    returnV = this->log_likelihood.bmd_hybrid_extra_bound(theta, BMD, BMRF,
                                                          isInc, advP);
    break;
  case CONTINUOUS_BMD_HYBRID_ADDED:
    break;
  default:
    break;
  }

  return returnV;
}

/*
 * Function: double returnBMD(Eigen::MatrixXd theta, contbmd BMDType, double
 * BMRF); Purpose: Depending on the type of BMD needed BMDType, the appropriate
 * BMD calculation algorithm is called given the parameters theta.
 *
 * Output: Returns the BMD of the given BMDtype. Returns a -1.0 BMD if the
 * BMDType is not defined
 */
template <class LL, class PR>
double cBMDModel<LL, PR>::returnBMD(Eigen::MatrixXd theta, contbmd BMDType,
                                    double BMRF, double advP) {

  // set all the fixed parameters to their value
  for (size_t i = 0; i < this->isFixed.size(); i++) {
    if (this->isFixed[i]) {
      theta(i, 0) = this->fixedV[i];
    }
  }

  double returnV = 0.0;

  switch (BMDType) {
  case CONTINUOUS_BMD_ABSOLUTE:
    returnV = this->log_likelihood.bmd_absolute(theta, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_STD_DEV:
    returnV = this->log_likelihood.bmd_stdev(theta, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_REL_DEV:
    returnV = this->log_likelihood.bmd_reldev(theta, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_POINT:
    returnV = this->log_likelihood.bmd_point(theta, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_EXTRA:
    returnV = this->log_likelihood.bmd_extra(theta, BMRF, isInc);
    break;
  case CONTINUOUS_BMD_HYBRID_EXTRA:
    returnV = this->log_likelihood.bmd_hybrid_extra(theta, BMRF, isInc, advP);
    break;
  case CONTINUOUS_BMD_HYBRID_ADDED:
    break;
  default:
    break;
  }

  return returnV;
}

/////////////////////////////////////////////////////////////////
// Function: equality_constraint(unsigned n,
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
// compute the negative penalized likelihood
//////////////////////////////////////////////////////////////////
template <class LL, class PR>
double cequality_constraint(unsigned n, const double *b, double *grad,
                            void *data) {

  c_optimInfo<LL, PR> *model = (c_optimInfo<LL, PR> *)data;
  Eigen::MatrixXd theta(n, 1);
  ////
  for (unsigned int i = 0; i < n; i++)
    theta(i, 0) = b[i]; // matrix copy
  ////

  // for the continuous models they are all equality constraints
  // thus we call equality and let the class cBMDstatmod figure out
  // what the equality constraint is
  return model->sm->equality_constraint(
      theta, grad, model->BMDType, model->cBMD, model->BMR, model->add_info);
}

/***************************************************************************************
cfindMAX_W_EQUALITY
Purpose: Does the necessary compuation to find the constraints
of the model given some set of equality constraints.
Input:
cBMDModel<LL, PR> *M - Pointer to a continous BMD model with likelihood LL, and
Prior PR Eigen::MatrixXd      - Starting values that will be fed into an
aditional starting value algorithm to start the optimizer. BMDType - The BMD
calculation type. BMD		     - Current BMD necessary for the equality
constraint BMRF                 - BMRF that is defined by the type of BMD (i.e.,
BMDType) isInc                - true if the curve is increasing, false
otherwise. tail_prob            - Used by the Hybrid approach only to define
"adverse" response.
*************************************************************************************/
// findMaxW_EQUALITY <- for continuous bmd models
template <class LL, class PR>
optimizationResult cfindMAX_W_EQUALITY(cBMDModel<LL, PR> *M,
                                       Eigen::MatrixXd start, contbmd BMDType,
                                       double BMD, double BMRF, bool isInc,
                                       double tail_prob) {
  // cout << start << endl << endl;
  int opt_iter = 0; // try several different optimization methods
                    // opt_iter starts of
  bool good_opt = false;
  // bool good_opt2 = true;
  optimizationResult oR;
  double minf;
  // double minf2;
  nlopt::result result;
  // nlopt::result result2;
  std::vector<double> x(start.rows());
  ///////////////////////////////////////////////////////////
  Eigen::MatrixXd temp_data = M->parmLB();
  std::vector<double> lb(M->nParms());
  for (int i = 0; i < M->nParms(); i++)
    lb[i] = temp_data(i, 0);
  temp_data = M->parmUB();
  std::vector<double> ub(M->nParms());
  for (int i = 0; i < M->nParms(); i++)
    ub[i] = temp_data(i, 0);

  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  start =
      M->starting_value(start, BMDType, BMD, BMRF, isInc, tail_prob, lb, ub);
  ///////////////////////////////////////////////////////////
  // cout << "equality start" << start << endl;

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = start(i, 0);
  }

  ///////////////////////////////////////////////////////////////////////////////
  c_optimInfo<LL, PR> info;
  info.BMR = BMRF;
  info.cBMD = BMD;
  info.sm = M;
  info.BMDType = BMDType;
  info.add_info = tail_prob;
  ///////////////////////////////////////////////////////////////////////////////

  while (opt_iter < 2 && !good_opt) {
    nlopt::opt opt(nlopt::LD_AUGLAG, M->nParms());      // alternate optimizer
    nlopt::opt local_opt(nlopt::LD_LBFGS, M->nParms()); // BOBYQA
    nlopt::opt local_opt2(nlopt::LN_SBPLX, M->nParms());
    /////////////////////////////////////////////////////////
    local_opt.set_xtol_abs(5e-5);
    local_opt2.set_xtol_abs(5e-5);
    local_opt.set_initial_step(5e-5);
    local_opt2.set_initial_step(5e-5);
    local_opt.set_maxeval(10000);
    local_opt2.set_maxeval(10000);
    /////////////////////////////////////////////////////////
    local_opt.set_lower_bounds(lb);
    local_opt.set_upper_bounds(ub);
    local_opt2.set_lower_bounds(lb);
    local_opt2.set_upper_bounds(ub);

    if (opt_iter == 0)
      opt.set_local_optimizer((const nlopt::opt)local_opt);
    else
      opt.set_local_optimizer((const nlopt::opt)local_opt2);
    /////////////////////////////////////////////////////////
    opt.add_equality_constraint(cequality_constraint<LL, PR>, &info, 1e-4);
    opt.set_min_objective(neg_pen_likelihood<LL, PR>, M);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_xtol_abs(5e-5);
    opt.set_maxeval(20000);
    ///////////////////////////////////////////////

    opt_iter++;              // iterate the optimization try counter
    result = nlopt::FAILURE; // Avoid uninit var exception checking result after
                             // NLOPT exception
    try {
      result = opt.optimize(x, minf);
      good_opt = true;
    } catch (nlopt::roundoff_limited &exc) {
      good_opt = false;
      // cout << "Error Round off" << endl;
    } catch (nlopt::forced_stop &exc) {
      good_opt = false;
      // cout << "Error Forced stop" << endl;
    } catch (const std::invalid_argument &exc) {
      good_opt = false;
      //	cout << "SHIT!!" << endl;
    } catch (const std::exception &exc) {
      good_opt = false;
      // cout << "Exception!!" << endl;
    }
    if (result > 5) { // Either 5 =
      // cout << "I am here:" <<  endl;
      good_opt = false;
    }
  }
  /*
  Eigen::Map<Eigen::MatrixXd> z(x.data(), M->nParms(), 1);
  Eigen::Map<Eigen::MatrixXd> zz(x2.data(), M->nParms(), 1);

  if (minf2 <= minf) {

          x = x2;
          minf = minf2;
          result = result2;
          good_opt = good_opt2;
  }
  */
  if (good_opt) { // if the opimization criteria worked
    Eigen::Map<Eigen::MatrixXd> d(x.data(), M->nParms(), 1); // return values
    oR.result = result;
    oR.functionV = minf;
    oR.max_parms = d;
  } else {
    oR.result = result;
    oR.functionV = std::numeric_limits<double>::quiet_NaN();
    oR.max_parms = Eigen::MatrixXd::Zero(M->nParms(), 1);
  }

  return oR;
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
// compute the negative penalized likelihood
//////////////////////////////////////////////////////////////////
template <class LL, class PR>
double neg_pen_likelihood_contbound(unsigned n, const double *b, double *grad,
                                    void *data) {

  c_optimInfo<LL, PR> *model = (c_optimInfo<LL, PR> *)data;
  model->sm->parameter_to_remove(model->BMDType);
  unsigned int p_remove = model->sm->parameter_to_remove(model->BMDType);
  int count = 0;
  //	double *tempgrad = new double[n + 1];
  std::vector<double> x(n + 1); // add one to the number of parameters
                                // LCO DEBUG ofstream file;
                                // file.open("bmds.log", fstream::app);
  // file << " ***HELLO" __FUNCTION__ << " at line: " << __LINE__ << endl;
  // file.close();
  // copy over the extra
  for (unsigned int i = 0; i < n + 1; i++) {
    if (i != p_remove) {
      x[i] = b[count];
      count++;
    }
  }
  count = 0;

  x = model->sm->bound_fix(x, model->BMDType, model->BMR, model->add_info,
                           model->cBMD, model->isInc);
  // cout << endl;
  // for (double b : x) cout << b << " ";
  // cout << endl;

  Eigen::MatrixXd theta(x.size(), 1);

  for (unsigned int i = 0; i < n + 1; i++) {
    theta(i, 0) = x[i];
  }
  // cout << model->sm->negPenLike(theta) << endl;
  if (grad) { // if there is a gradient
    Eigen::MatrixXd mgrad = model->sm->gradient(theta);

    for (unsigned int i = 0; i < n + 1; i++) {
      if (i != p_remove) {
        grad[count] = mgrad(i, 0); // remove the
        count++;
      }
    }
  }
  //	delete tempgrad;

  return model->sm->negPenLike(theta);
}

/***************************************************************************************
cfindMAX_W_BOUND
Purpose: Does the necessary compuation to find the constraints
of the model given that the BMD is computed based upon the bound of the BMD
Input:
cBMDModel<LL, PR> *M - Pointer to a continous BMD model with likelihood LL, and
Prior PR Eigen::MatrixXd      - Starting values that will be fed into an
aditional starting value algorithm to start the optimizer. BMDType - The BMD
calculation type. BMD		     - Current BMD necessary for the equality
constraint BMRF                 - BMRF that is defined by the type of BMD (i.e.,
BMDType) isInc                - true if the curve is increasing, false
otherwise. tail_prob            - Used by the Hybrid approach only to define
"adverse" response.
*************************************************************************************/
// findMaxW_EQUALITY <- for continuous bmd models
template <class LL, class PR>
optimizationResult cfindMAX_W_BOUND(cBMDModel<LL, PR> *M, Eigen::MatrixXd start,
                                    contbmd BMDType, double BMD, double BMRF,
                                    bool isInc, double tail_prob) {
  optimizationResult oR;
  double minf = 0;
  nlopt::result result = nlopt::FAILURE;
  int vecSize = start.rows() - 1;
  std::vector<double> x(vecSize); // drop the number of parameters by 1
  std::vector<double> lb(vecSize);
  std::vector<double> ub(vecSize);
  Eigen::MatrixXd datal = M->parmLB();
  Eigen::MatrixXd datau = M->parmUB();
  int count = 0;
  ///////////////////////////////////////////////////////////////////////////////
  // remove the extra parameter from the list
  ///////////////////////////////////////////////////////////////////////////////
  int p_remove = M->parameter_to_remove(BMDType);
  for (int i = 0; i < M->nParms(); i++) {
    if (i != p_remove) {
      lb[count] = datal(i, 0);
      ub[count] = datau(i, 0);
      double temp = start(i, 0);
      if (temp < lb[count])
        temp = lb[count];
      else if (temp > ub[count])
        temp = ub[count];
      x[count] = temp;
      count++;
    } // end if
  }   // end for

  ///////////////////////////////////////////////////////////////////////////////
  c_optimInfo<LL, PR> info;
  info.BMR = BMRF;
  info.cBMD = BMD;
  info.isInc = isInc;
  info.sm = M;
  info.BMDType = BMDType;
  info.add_info = tail_prob;

  bool good_opt = false;
  int opt_iter = 0;
  //////////////////////////////////////////////////////////////////////////////
  nlopt::opt opt2(nlopt::LD_LBFGS, vecSize);
  opt2.set_initial_step(1e-4);
  opt2.set_min_objective(neg_pen_likelihood_contbound<LL, PR>, &info);
  opt2.set_lower_bounds(lb);
  opt2.set_upper_bounds(ub);
  opt2.set_xtol_abs(5e-4);
  opt2.set_maxeval(20000);
  //////////////////////////////////////////////////////////////////////////////
  nlopt::opt opt(nlopt::LN_BOBYQA, vecSize);
  opt.set_initial_step(1e-4);
  opt.set_min_objective(neg_pen_likelihood_contbound<LL, PR>, &info);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_abs(5e-4);
  opt.set_maxeval(20000);

  nlopt::opt opt3(nlopt::LN_SBPLX, vecSize);
  opt3.set_initial_step(1e-4);
  opt3.set_min_objective(neg_pen_likelihood_contbound<LL, PR>, &info);
  opt3.set_lower_bounds(lb);
  opt3.set_upper_bounds(ub);
  opt3.set_xtol_abs(5e-4);
  opt3.set_maxeval(20000);
  ///////////////////////////////////////////////////////////////////////////////

  //	if(M->modelling_type() == cont_model::gamma_aerts){
  //		opt_iter = 1;
  //	}
  //    cout << " cFINDMAX init=" << endl;
  //    for (auto i: x)
  //        std::cout << i << ' ';

  while (!good_opt && opt_iter <= 2) {
    try {

      switch (opt_iter) {
      case 0:
        opt_iter++;
        result =
            opt2.optimize(x, minf); // fastest algorithm first
                                    //				cout << "huh1";
        break;
      case 1:
        opt_iter++;
        result =
            opt3.optimize(x, minf); // second fastest algorithm next
                                    //				cout << "huh2";
        break;
      default:
        opt_iter++;
        result = opt.optimize(
            x,
            minf); // most stable one third -- but slower
                   //				cout << "huh3 " << result;
      }
      // file << "result= " << result << ", minf= " << minf << endl;
      // flush(file);
      good_opt = true;
      // opt_iter++;
    } catch (nlopt::roundoff_limited &exc) {
      good_opt = false;
      DEBUG_LOG(file, "opt_iter= " << opt_iter << ", error: roundoff_limited");
      //	cout << "Error Round off" << endl;
    } catch (nlopt::forced_stop &exc) {
      good_opt = false;
      DEBUG_LOG(file, "opt_iter= " << opt_iter << ", error: roundoff_limited");
      //	cout << "Error Forced stop" << endl;
    } catch (const std::invalid_argument &exc) {
      good_opt = false;
      DEBUG_LOG(file, "opt_iter= " << opt_iter
                                   << ", error: invalid arg: " << exc.what());
    } catch (const std::exception &exc) {
      good_opt = false;
      DEBUG_LOG(file,
                "opt_iter= " << opt_iter << ", general error: " << exc.what());
      // cout << "Exception!!" << endl;
    }
    if (result >= 5) { // Either 5 =
      good_opt = false;
    }
    DEBUG_LOG(file, "\topt_iter= " << opt_iter << ", result= " << result
                                   << ", minf= " << minf
                                   << ", good_opt= " << good_opt);
  }
  //	cout << "Opt "<< good_opt << endl;
  Eigen::MatrixXd xxx = Eigen::MatrixXd::Zero(M->nParms(), 1);
  // std::vector<double> xxx(x.size() + 1);
  count = 0;
  for (int i = 0; i < M->nParms(); i++) {
    if (i != p_remove) {
      xxx(i, 0) = x[count];
      count++;
    }
  }

  if (good_opt) { // if the opimization criteria worked
    x = M->bound_fix(xxx, BMDType, BMRF, tail_prob, BMD, isInc);
    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(x.size(), 1);
    for (int i = 0; i < x.size(); ++i) {
      result_matrix(i, 0) = x[i]; 
    }
    oR.result = result;
    oR.functionV = minf;
    oR.max_parms = result_matrix;
  } else {
    oR.result = result;
    oR.functionV = std::numeric_limits<double>::quiet_NaN();
    oR.max_parms = Eigen::MatrixXd::Zero(M->nParms(), 1);
  }

  DEBUG_CLOSE_LOG(file);

  return oR;
}

///////////////////////////////////////////////////////////////////////////////
// Function profile_BMDNC(dBMDModel<LL, PR>  *M,
//						 bool isExtra,		// true
// if it is
// false if it is added 						 double
// BMR, double BMDchange, double totalChange, bool robust) Purpose: This
// function iteratively changes the BMD by a BMDchange%
//          until a total change in the penalized likelihood is found.
// Input  : dBMDModel<LL, PR>  *M - Dichotomous BMD model
//		   bool isExtra          - true if it is extra risk
//          double BMD			 - current BMD at the MAP
//		   double BMR            - BMR [0,1]
//	       double BMDchange      - %Change in the BMD to increment each time
//		   double totalChange    - totalChange in penalized likelihood
// before one stops
//          bool   robust         - true if we do a robust search of the
//          optimization space, false otherwise
////////////////////////////////////////////////////////////////////////////////
template <class LL, class PR>
Eigen::MatrixXd profile_cBMDNC(cBMDModel<LL, PR> *M, contbmd BMDType,
                               const double BMD, const double BMRF,
                               const double tail_prob, const double BMDchange,
                               const double totalChange, bool isIncreasing) {

  // double mapBMD = BMD; // current BMD evaluated at estimate M->getEST
  Eigen::MatrixXd parms = M->getEST();
  Eigen::MatrixXd result(3, 1);
  optimizationResult oR;

  double MAP_LIKE = M->negPenLike(parms);
  Eigen::MatrixXd ret1(3, 1), ret2, fParms;
  std::list<Eigen::MatrixXd> CL, temp1, temp2;

  double CLIKE = MAP_LIKE;
  int max_iters = 0;
  double CBMD = BMD * (1.0 - BMDchange);
  double PLIKE = CLIKE;
  // double PBMD = BMD;
  bool error = false;

  ret1(0, 0) = PLIKE;
  ret1(1, 0) = BMD;
  ret1(2, 0) = 666;
  CL.push_front(ret1);

  //	cout << "BMDchange : " << BMDchange << " : " <<
  // M->type_of_profile(BMDType) << endl;

  while (fabs(MAP_LIKE - CLIKE) < totalChange && max_iters < 300) {

    // fit the profile likelihood
    if (M->type_of_profile(BMDType) == PROFILE_EQUALITY) {
      oR = cfindMAX_W_EQUALITY<LL, PR>(M, parms, BMDType, CBMD, BMRF,
                                       isIncreasing, tail_prob);
    } else {
      oR = cfindMAX_W_BOUND<LL, PR>(M, parms, BMDType, CBMD, BMRF, isIncreasing,
                                    tail_prob);
    }

    parms = oR.max_parms;
    result(0, 0) = oR.functionV;
    result(1, 0) = CBMD;
    result(2, 0) = oR.result;
    temp1.push_front(parms);

    CLIKE = oR.functionV;
    PLIKE = CLIKE;
    CBMD *= 1 - BMDchange;
    CL.push_front(result);
    max_iters++;

    if (std::isnan(oR.functionV)) {
      error = true;
      //			cout << "nans 1 ";
    }
  }

  CLIKE = MAP_LIKE;
  CBMD = BMD * (1.0 + BMDchange);
  max_iters = 0;
  parms = M->getEST();
  error = false;

  while (fabs(MAP_LIKE - CLIKE) < totalChange && max_iters < 300 && !error) {
    if (M->type_of_profile(BMDType) == PROFILE_EQUALITY) {
      oR = cfindMAX_W_EQUALITY<LL, PR>(M, parms, BMDType, CBMD, BMRF,
                                       isIncreasing, tail_prob);
    } else {
      oR = cfindMAX_W_BOUND<LL, PR>(M, parms, BMDType, CBMD, BMRF, isIncreasing,
                                    tail_prob);
    }

    parms = oR.max_parms;
    result(0, 0) = oR.functionV;
    result(1, 0) = CBMD;
    result(2, 0) = oR.result;
    CLIKE = oR.functionV;
    PLIKE = CLIKE;
    CBMD *= 1 + BMDchange;

    if (std::isnan(oR.functionV) || std::isinf(CBMD)) {
      error = true;
      //			cout << std::isnan(oR.functionV) <<
      // std::isinf(CBMD) << "nan/inf2 ";
    }

    CL.push_back(result);

    max_iters++;
  }

  // Take the lists and put them in one big
  // matrix to return.
  Eigen::MatrixXd returnMat(CL.size(), 3);
  int ii = 0;
  for (Eigen::MatrixXd it : CL) {
    returnMat.row(ii) = it.transpose();
    ii++;
  }

  returnMat.col(0).array() =
      round(round(1e4 * returnMat.col(0).array()) - round(1e4 * MAP_LIKE)) /
      1e4; // round it to be within the optimizers tol

  //	cout << returnMat << endl;

  return returnMat;
}

#endif
