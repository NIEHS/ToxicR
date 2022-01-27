/* Copyright 2021  NIEHS <matt.wheeler@nih.gov>
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
#include "bmd_calculate.h"
#include "bmdStruct.h"

#if defined R_COMPILATION || defined __cplusplus
  #include <cmath> 
#else
  #include <math.h>
#endif

#ifdef R_COMPILATION  
  #include <RcppEigen.h>
  #include <RcppGSL.h>
  using namespace Rcpp;
  using Rcpp::as;
#elif __cplusplus
  #include <Eigen/Dense>
#endif
  
#define NUM_PRIOR_COLS 5

#include <algorithm>
#include <vector>
#include <limits>
#include <set>


#include "lognormalTests.h"
#include "normalTests.h"
#include <iostream>
/* log_normal_AOD
 * 
 * 
 * 
 * 
 */
void log_normal_AOD_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, 
                         bool bSuffStat, continuous_deviance * CD){

  /////////////////////////////////////////////////////////////////////////////////////
  // Test A1
  lognormalLLTESTA1 a1Test(Y, X, true);
  // if it is a sufficient statistics model
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  std::vector<double> udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd meanX = Eigen::MatrixXd::Zero(Y.rows(), udoses.size()); 
  Eigen::MatrixXd meanX2 = Eigen::MatrixXd::Ones(Y.rows(),1);
  Eigen::MatrixXd W     = Eigen::MatrixXd::Zero(Y.rows(), Y.rows()); 
  for (int i = 0; i < Y.rows(); i++){
    W(i,i) = Y(i,2); 
  }
  for (int i = 0; i < meanX.rows(); i++)
  {
    for (int j = 0; j < udoses.size(); j++) {
      meanX(i, j) = udoses[j] == X(i, 0) ? 1.0 : 0.0; 
    }
  }
  Eigen::MatrixXd InitA = (meanX.transpose()*W*meanX);
  Eigen::MatrixXd InitB = (meanX2.transpose()*W*meanX2);
  InitA = InitA.inverse()*meanX.transpose()*W*Y;
  InitB = InitB.inverse()*meanX2.transpose()*W*Y;
  //Eigen::MatrixXd startV1(a1Test.nParms(),1); 
  
  int nParms = a1Test.nParms();
  std::vector<double> fix1(nParms); for (unsigned int i = 0; i < nParms; i++) { fix1[i] = 0.0; }
  std::vector<bool> isfix1(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix1[i] = false; }
  Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a1Priors.row(i) << 0, 0, 1, 0, 1e7;
  }
  a1Priors.row(nParms - 1) << 0, 0, 1, -1e7, 1e7;
  
  Eigen::MatrixXd startV1(a1Test.nParms(),1); 
  double tempV = 0.0; 
  for (unsigned int i = 0; i < InitA.rows() - 1; i++){
    startV1(i,0) = exp(InitA(i,0));
  }
 
  startV1(a1Test.nParms() - 1,0) = InitB(0,1); 
   
  IDcontinuousPrior a1Init(a1Priors);
  statModel<lognormalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                          isfix1, fix1);
  optimizationResult a1Result = findMAP<lognormalLLTESTA1, IDcontinuousPrior>(&a1Model,startV1,0);
  
  CD->A1 = a1Result.functionV;
  CD->N1 = a1Result.max_parms.rows();
  /////////////////////////////////////////////////////////////////////////////////////
  // Test A2
  lognormalLLTESTA2 a2Test(Y, X, bSuffStat);
  nParms = a2Test.nParms();
  std::vector<double> fix2(nParms); for (unsigned int i = 0; i < nParms; i++) { fix2[i] = 0.0; }
  std::vector<bool> isfix2(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix2[i] = false; }
  Eigen::MatrixXd a2Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  for (unsigned int i = 0; i < nParms / 2; i++) {
    a2Priors.row(i) << 0, a1Result.max_parms(i), 1, 0, 1e8;
  }
  IDcontinuousPrior a2Init(a2Priors);
  Eigen::MatrixXd startV2(nParms,1);
  for (unsigned int i = 0; i < nParms; i++) {
    a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  for (unsigned int i = 0; i < nParms / 2; i++) {
    a2Priors.row(i) << 0, a1Result.max_parms(i), 1, -1e8, 1e8;
  }
  
  statModel<lognormalLLTESTA2, IDcontinuousPrior> a2Model(a2Test, a2Init,
                                                          isfix2, fix2);
  optimizationResult a2Result = findMAP<lognormalLLTESTA2, IDcontinuousPrior>(&a2Model,startV2,0);
  CD->A2 = a2Result.functionV;
  CD->N2 = a2Result.max_parms.rows();
  /////////////////////////////////////////////////////////////////////////////////////
  // Test R
  lognormalLLTESTR rTest(Y, X, bSuffStat);
  nParms = rTest.nParms();
  std::vector<double> fix4(nParms); for (unsigned int i = 0; i < nParms; i++) { fix4[i] = 0.0; }
  std::vector<bool> isfix4(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix4[i] = false; }
  Eigen::MatrixXd rPriors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    rPriors.row(i) << 0, 1, 1, 0, 1e8;
  }
  rPriors.row(nParms - 1) << 0, 0, 1, -1e8, 1e8;

  Eigen::MatrixXd startR(2,1);
  startR(0,0) = exp(InitB(0,0));
  startR(1,0) = InitB(0,1);
  IDcontinuousPrior rInit(rPriors);
  statModel<lognormalLLTESTR, IDcontinuousPrior> rModel(rTest, rInit,
                                                        isfix4, fix4);
  optimizationResult rResult = findMAP<lognormalLLTESTR, IDcontinuousPrior>(&rModel,startR,0);
  
  CD->R = rResult.functionV;
  CD->NR = rResult.max_parms.rows();
}

/* normal_AOD
 * 
 * 
 * 
 */
void normal_AOD_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, 
                     bool bSuffStat, continuous_deviance * CD){
  
  //////////////////////////////////////////////////////////////////////

  normalLLTESTA1 a1Test(Y, X, true);
  int  nParms = a1Test.nParms();
  std::vector<double> fix1(nParms); for (unsigned int i = 0; i < nParms; i++) { fix1[i] = 0.0; }
  std::vector<bool> isfix1(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix1[i] = false; }
  Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a1Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  IDcontinuousPrior a1Init(a1Priors);
  
  // if it is a sufficient statistics model
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  std::vector<double> udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd meanX = Eigen::MatrixXd::Zero(Y.rows(), udoses.size()); 
  Eigen::MatrixXd meanX2 = Eigen::MatrixXd::Ones(Y.rows(),1);
  Eigen::MatrixXd W     = Eigen::MatrixXd::Zero(Y.rows(), Y.rows()); 
 
  for (int i = 0; i < Y.rows(); i++){
    W(i,i) = Y(i,2); 
  }
 
  for (int i = 0; i < meanX.rows(); i++)
  {
    for (int j = 0; j < udoses.size(); j++) {
      meanX(i, j) = udoses[j] == X(i, 0) ? 1.0 : 0.0; 
    }
  }
 
  Eigen::MatrixXd InitA = (meanX.transpose()*W*meanX);
  Eigen::MatrixXd InitB = (meanX2.transpose()*W*meanX2);
  InitA = InitA.inverse()*meanX.transpose()*W*Y;
  InitB = InitB.inverse()*meanX2.transpose()*W*Y;
 
  Eigen::MatrixXd startV1(a1Test.nParms(),1); 
  
  double tempV = 0.0; 
 
  for (unsigned int i = 0; i < a1Test.nParms()-1; i++){
    startV1(i,0) = InitA(i,0);
  }
  startV1(a1Test.nParms() - 1,0) = log(InitB(0,1)*InitB(0,1)); 
 
  statModel<normalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                       isfix1, fix1);
 
  optimizationResult a1Result = findMAP<normalLLTESTA1, IDcontinuousPrior>(&a1Model,startV1,0);
  
  CD->A1 = a1Result.functionV;
  CD->N1 = a1Result.max_parms.rows(); 
  ////////////////////////////////////////////////////////////////////
  normalLLTESTA2 a2Test(Y, X, bSuffStat);
  nParms = a2Test.nParms();
  std::vector<double> fix2(nParms); for (unsigned int i = 0; i < nParms; i++) { fix2[i] = 0.0; }
  std::vector<bool> isfix2(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix2[i] = false; }
  Eigen::MatrixXd a2Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  for (unsigned int i = 0; i < nParms / 2; i++) {
    a2Priors.row(i) << 0, a1Result.max_parms(i), 1, -1e8, 1e8;
  }
  
  IDcontinuousPrior a2Init(a2Priors);
  Eigen::MatrixXd startV2(nParms,1);
  for (unsigned int i = 0; i < nParms/2; i++ ){
    startV2(i,0) = InitA(i,0); 
  } 
  
  for (unsigned int i = 0;  i < nParms/2; i++ ){
    startV2(i+nParms/2,0) = log(InitA(i,1)*InitA(i,1)); 
  }
 
  statModel<normalLLTESTA2, IDcontinuousPrior> a2Model(a2Test, a2Init,
                                                       isfix2, fix2);

  optimizationResult a2Result = findMAP<normalLLTESTA2, IDcontinuousPrior>(&a2Model,startV2,0);
  CD->A2 = a2Result.functionV;
  CD->N2 = a2Result.max_parms.rows(); 
  //////////////////////////////////////////////////////////////////////
  normalLLTESTA3 a3Test(Y, X, bSuffStat);
  nParms = a3Test.nParms();
  std::vector<double> fix3(nParms); for (unsigned int i = 0; i < nParms; i++) { fix3[i] = 0.0; }
  std::vector<bool> isfix3(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix3[i] = false; }
  Eigen::MatrixXd a3Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a3Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
 

  a3Priors.row(nParms - 2) << 0, a2Result.max_parms(nParms - 2), 1, -1e4, 1e4;
  a3Priors.row(nParms - 1) << 0, a2Result.max_parms(nParms - 1), 0, -1e4, 1e4;
 
  IDcontinuousPrior a3Init(a3Priors);

  Eigen::MatrixXd startV3(nParms,1);
  Eigen::MatrixXd meanX3 = Eigen::MatrixXd::Ones(InitA.rows(),2);
  for (int i = 0; i < meanX3.rows(); i++){
    meanX3(i,1) = InitA(i,0) > 0? log(InitA(i,0)):-13.8;  
    InitA(i,1) = log(InitA(i,1)*InitA(i,1));   
  }
 
  for (unsigned int i = 0; i < nParms-2; i++ ){
    startV3(i,0) = a2Result.max_parms(i); 
  } 


  Eigen::MatrixXd InitC = (meanX3.transpose()*meanX3); 
  InitC = InitC.inverse()*meanX3.transpose()*InitA.col(1); 
  startV3(nParms-2,0) = 0;
  startV3(nParms-1,0) = a1Result.max_parms(a1Result.max_parms.rows()-1,0);

  
  statModel<normalLLTESTA3, IDcontinuousPrior> a3Model(a3Test, a3Init,
                                                       isfix3, fix3);
  optimizationResult a3Result = findMAP<normalLLTESTA3, IDcontinuousPrior>(&a3Model,startV3,0);
  CD->A3 = a3Result.functionV;
  CD->N3 = a3Result.max_parms.rows(); 
  
  //////////////////////////////////////////////////////////////
  normalLLTESTR RTest(Y, X, bSuffStat);
  nParms = RTest.nParms();
  std::vector<double> fixR(nParms); for (unsigned int i = 0; i < nParms; i++) { fixR[i] = 0.0; }
  std::vector<bool> isfixR(nParms); for (unsigned int i = 0; i < nParms; i++) { isfixR[i] = false; }
  Eigen::MatrixXd RPriors(nParms, NUM_PRIOR_COLS);

  for (unsigned int i = 0; i < nParms; i++) {
     RPriors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  IDcontinuousPrior RInit(RPriors);
  Eigen::MatrixXd startR(2,1);
  startR(0,0) = InitB(0,0);
  startR(1,0) = log(InitB(0,1)*InitB(0,1));
  statModel<normalLLTESTR, IDcontinuousPrior> RModel(RTest, RInit,
                                                       isfixR, fixR);
  optimizationResult RResult = findMAP<normalLLTESTR, IDcontinuousPrior>(&RModel,startR);
  CD->R = RResult.functionV;
  CD->NR = RResult.max_parms.rows(); 
}

void variance_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, bool bSuffStat,
                          double *v_c, double *v_nc, double *v_pow){
  //////////////////////////////////////////////////////////////////////
  
  normalLLTESTA1 a1Test(Y, X, true);
  int  nParms = a1Test.nParms();
  std::vector<double> fix1(nParms); for (unsigned int i = 0; i < nParms; i++) { fix1[i] = 0.0; }
  std::vector<bool> isfix1(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix1[i] = false; }
  Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
 
  for (unsigned int i = 0; i < nParms; i++) {
    a1Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  IDcontinuousPrior a1Init(a1Priors);
  
  // if it is a sufficient statistics model
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  std::vector<double> udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd meanX = Eigen::MatrixXd::Zero(Y.rows(), udoses.size()); 
  Eigen::MatrixXd meanX2 = Eigen::MatrixXd::Ones(Y.rows(),1);
  Eigen::MatrixXd W     = Eigen::MatrixXd::Zero(Y.rows(), Y.rows()); 

  for (int i = 0; i < Y.rows(); i++){
    W(i,i) = Y(i,2); 
  }

  for (int i = 0; i < meanX.rows(); i++)
  {
    for (int j = 0; j < udoses.size(); j++) {
      meanX(i, j) = udoses[j] == X(i, 0) ? 1.0 : 0.0; 
    }
  }
  Eigen::MatrixXd InitA = (meanX.transpose()*W*meanX);
  Eigen::MatrixXd InitB = (meanX2.transpose()*W*meanX2);
  InitA = InitA.inverse()*meanX.transpose()*W*Y;
  InitB = InitB.inverse()*meanX2.transpose()*W*Y;
  
  Eigen::MatrixXd startV1(a1Test.nParms(),1); 
  
  double tempV = 0.0; 
 
  for (unsigned int i = 0; i < a1Test.nParms() - 1; i++){
    startV1(i,0) = InitA(i,0);
  }
  startV1(a1Test.nParms() - 1,0) = InitB(0,1)*InitB(0,1); 

  statModel<normalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                       isfix1, fix1);
  
  optimizationResult a1Result = findMAP<normalLLTESTA1, IDcontinuousPrior>(&a1Model,startV1,0);
  
  *v_c = a1Result.max_parms(a1Result.max_parms.rows()-1,0); 
  //////////////////////////////////////////////////////////////////////////////////////////////

  normalLLTESTA3 a3Test(Y, X, bSuffStat);
  nParms = a3Test.nParms();
  std::vector<double> fix3(nParms); for (unsigned int i = 0; i < nParms; i++) { fix3[i] = 0.0; }
  std::vector<bool> isfix3(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix3[i] = false; }
  Eigen::MatrixXd a3Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a3Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  
  for (unsigned int i = 0; i < nParms - 2; i++) {
    a3Priors.row(i) << 0, 0, 1, 0, 1e8;
  }
  IDcontinuousPrior a3Init(a3Priors);
  
  Eigen::MatrixXd startV3(nParms,1);
  Eigen::MatrixXd meanX3 = Eigen::MatrixXd::Ones(InitA.rows(),2);
  for (int i = 0; i < meanX3.rows(); i++){
    meanX3(i,1) = InitA(i,0) > 0? log(InitA(i,0)):-13.8;  
    InitA(i,1) = log(InitA(i,1)*InitA(i,1));   
  }
  for (unsigned int i = 0; i < nParms-2; i++ ){
    startV3(i,0) = InitA(i,0); 
  } 

  Eigen::MatrixXd InitC = (meanX3.transpose()*meanX3); 
  InitC = InitC.inverse()*meanX3.transpose()*InitA.col(1); 
  startV3(nParms-2) = InitC(1,0);
  startV3(nParms-1) = InitC(0,0);
  
  statModel<normalLLTESTA3, IDcontinuousPrior> a3Model(a3Test, a3Init,
                                                       isfix3, fix3);
  optimizationResult a3Result = findMAP<normalLLTESTA3, IDcontinuousPrior>(&a3Model,startV3);
  *v_nc  = a3Result.max_parms(a3Result.max_parms.rows()-1,0);
  *v_pow =   a3Result.max_parms(a3Result.max_parms.rows()-2,0);

 
  return; 
}
