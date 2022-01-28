#######################################################################################################
# Copyright 2020  NIEHS <matt.wheeler@nih.gov>
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
#is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################################################
.fit_types                <- c('laplace','mle','mcmc')

#continuous global information
.continuous_models        <- c("hill","exp-3","exp-5","power","FUNL","polynomial")
.continuous_distributions <- c("normal","normal-ncv","lognormal")
.continuous_risk_types    <- c('abs','sd','rel','hybrid')

#dichotomous global information

.dichotomous_models       <- c("hill","gamma","logistic","log-logistic","log-probit","multistage","probit",
                               "qlinear","weibull")
