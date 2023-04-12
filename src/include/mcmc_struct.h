/*
 * Copyright 2020  US HHS, NIEHS 
 * Author Matt Wheeler 
 * e-mail: <matt.wheeler@nih.gov> 
 *
 *
 *Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 *and associated documentation files (the "Software"), to deal in the Software without restriction,
 *including without limitation the rights to use, copy, modify, merge, publish, distribute,
 *sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 *is furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all copies
 *or substantial portions of the Software.

 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 *CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 *OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 */
#ifndef _DEFINE_MCMC_STRUCT_H
#define _DEFINE_MCMC_STRUCT_H

#ifdef R_COMPILATION
// necessary things to run in R
#include <RcppEigen.h>
#include <RcppGSL.h>
#else
#include <Eigen/Dense>
#endif

class mcmcSamples
{
public:
  /*******************
   *for model averaging
   *********************/
  double map;
  Eigen::MatrixXd map_estimate;
  Eigen::MatrixXd map_cov;
  /**********************/
  Eigen::MatrixXd BMD;
  Eigen::MatrixXd samples;
  Eigen::MatrixXd log_posterior;
  double BMR;
  bool isExtra;
  mcmcSamples()
  {
    isExtra = true;
    map = 0;
    BMR = 0;
  }
};

#endif
