################################################
# Copyright 2021  NIEHS <matt.wheeler@nih.gov>
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
# 
#   
#################################################
.crutial_stat_constant <- function(param,y,doses,var,mean_function,decreasing,alpha=0){
        expected  <- mean_function(param,doses,decreasing)
         sq_resid <- (y-expected)^2/(var*expected^alpha)
        return(sum(sq_resid))
}

#@ .pvalue_cont_mcmc
#@ Function that computes a p-value from an MCMC object based upon 
#@ the method of Johnson(2007) on Bayesian pivitol quantities
.pvalue_cont_mcmc <- function(mcmc_fit,model,y,doses,distribution,decreasing){
  pValue_return = NA; 
 
    if (model == "FUNL"){
      func = .cont_FUNL_f
    }
    if (model == "exp-5"){
      func = .cont_exp_5_f
    }
    if (model == "exp-3"){
      func = .cont_exp_3_f
    }
    if (model == "hill"){
      func = .cont_hill_f
    }
    if (model == "power"){
      func = .cont_power_f
    }
    if (model == "polynomial"){
      
    }
  if (distribution == "normal"){
    q<- apply(mcmc_fit$mcmc_result$PARM_samples[1000:nrow(mcmc_fit$mcmc_result$PARM_samples),], 1,.crutial_stat_constant,y=y,
              doses=doses,var =exp(mcmc_fit$varOpt[1]),mean_function=func,decreasing=decreasing,alpha=0)
    temp <- pchisq(quantile(q,0.90),length(y)-1)
    pValue_return  = 1 - max(0,(temp*length(q)-0.90*length(q)+1)/(length(q)-0.90*length(q)+1))
  }
  
  if (distribution == "normal-ncv"){
    q<- apply(mcmc_fit$mcmc_result$PARM_samples[1000:nrow(mcmc_fit$mcmc_result$PARM_samples),], 1,.crutial_stat_constant,y=y,
              doses=doses,var =exp(mcmc_fit$varOpt[2]),mean_function=func,decreasing=decreasing,alpha=mcmc_fit$varOpt[3])
    temp <- pchisq(quantile(q,0.90),length(y)-1)
    pValue_return  = 1 - max(0,(temp*length(q)-0.90*length(q)+1)/(length(q)-0.90*length(q)+1))
  }
  
  return(pValue_return)
}