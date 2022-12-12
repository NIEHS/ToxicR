/*
* Copyright 2020  US. Department of Health and Human Services (HHS), 
 * National Institute of Environmental Health Sciences (NIEHS)
 * Email: Matt Wheeler  <matt.wheeler@nih.gov>
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
#ifndef __ANALYSIS_OF_DEVIANCE
#define __ANALYSIS_OF_DEVIANCE

void log_normal_AOD_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, 
                         bool bSuffStat, continuous_deviance * CD);

void normal_AOD_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, 
                     bool bSuffStat, continuous_deviance * CD);

void variance_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, bool bSuffStat,
                   double *v_c, double *v_nc, double *v_pow);

#endif



