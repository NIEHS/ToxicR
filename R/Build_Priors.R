#Copyright 2020  NIEHS <matt.wheeler@nih.gov>
#   
#
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

.parse_prior<-function(prior){
  rV <-list()
  rV$prior <- prior$prior
  
  temp_a  <- regexpr("[[][a-zA-Z]+-*[a-zA-Z]+[]]",prior$model)
  start   <- temp_a[1] + 1
  end     <- start + attr(temp_a,"match.length") - 3
  if(temp_a == -1){
    stop("Could not find a distribution for analysis.")
  }
  rV$distribution = substr(prior$model,start,end)
  rV$model = prior$mean
  return(rV)
  
}

#' @title create_continuous_prior  Given priorlist, a model, 
#'        and a distribution. Create a prior for a given analysis. 
#' @param prior_list First Prior
#' @param model Model to be used
#' one of \code{"hill","exp-3","exp-5","power","polynomial"}
#' @param distribution - Normal "normal", Normal non-constant variance "normal-ncv", or
#'                       log-normal "lognormal"
#' @param deg - For polynomial models only, the degree of the polynomial. 
#' @return new BMDprior list. This object is essentially a prior list constructed by 
#' \code{create_prior_lists} with a model type and variance type. 
#' 
#' @examples 
#' plist<- create_prior_list(normprior(0,0.1,-100,100), # a
#'                           normprior(0,1,  -1e2,1e2),     # b
#'                           lnormprior(1,0.2,0,18),  #k
#'                           normprior(0,1,-18,18))
#'  
#'  power_normal <- create_continuous_prior(plist,"power","normal") 
#'
create_continuous_prior <- function( prior_list,model,distribution,deg=2){

  if (!("BMDmodelprior" %in% class(prior_list))){
    stop("Prior is not of a 'BMDmodelprior' class. A probable solution is to 
          define the prior using function `create_prior_list`.")
  }
  if (!(model %in% .continuous_models )){
    stop(cat("Model Type must be one of:",.continuous_models,"\n"))
  }
  if (!(distribution %in% .continuous_distributions )){
    stop(cat("Distribution must be one of:",.continuous_distributions,"\n"))
  }
  temp = floor(deg)
  
  if ( deg < 1){
    stop("Polynomial degree must be greater than or equal to 1.")
  }
  if ( temp != deg){
    stop("Polynomial degree must be an integer.")
  }
  
  p = NA
  
  if ("hill" == model){
    p = .check_hill(prior_list,distribution) 
  }
  
  if ("FUNL" == model){
    p = .check_FUNL(prior_list,distribution) 
  }
 
  if ("exp-5" == model){
    p = .check_exp5(prior_list,distribution) 
  }
  
  if ("exp-3" == model){
    p = .check_exp3(prior_list,distribution) 
  }
  
  if ("polynomial" == model){
    p = .check_polynomial(prior_list,distribution) 
  }
  
  if ("power" == model){
    p = .check_power(prior_list,distribution) 
  }
  class(p)<- "BMD_Bayes_continuous_model"
 
  return(p)
}

.check_hill <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  if (distribution == "normal"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Normal Hill model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Hill Model [normal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) != 6){
      stop("Normal-NCV Hill model prior requires 6 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5,4] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Hill Model [normal-ncv]"
    prior$parameters <- c("a","b","c","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/Hill specification is not presently available.")
  }
  prior$mean = .continuous_models[1]
  return(prior)
}

.check_exp5 <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
  if (distribution == "normal"){
   
    if (nrow(temp) != 5){
      stop("Normal Exponential-5  model prior requires 5 parameters.")
    }
    if (sum(temp[,4] > temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    
    prior$model = "Exponential-5 Model [normal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    if (nrow(temp) != 6){
      stop("Normal Exponential-5 model prior requires 6 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5,4] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Exponential-5 [normal-ncv]"
    prior$parameters <- c("a","b","c","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Lognormal Exponential-5  model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Exponential-5 Model [lognormal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  prior$mean = .continuous_models[3]
  return(prior)
}

.check_power <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  if (distribution == "normal"){
    temp = prior[[1]]
    if (nrow(temp) != 4){
      stop("Normal Power model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[3,4] < 0){
      stop("The power parameter d (parameter 3) can not have a lower bound less than zero.")
    }
    prior$model = "Power Model [normal]"
    prior$parameters <- c("a","b","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Normal-NCV Power model prior requires 5 parameters.")
    }
    if (sum(temp[,4]> temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[3,4] < 0){
      stop("The power parameter d (parameter 3) can not have a lower bound less than zero.")
    }
    
    if (temp[4,4] < 0){
      stop("The prior on \rho (parameter 4) can not have a lower bound less than zero.")
    }
    
    prior$model = "Power Model [normal-ncv]"
    prior$parameters <- c("a","b","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/Power specification is not presently available.")
  }
  prior$mean = .continuous_models[4]
  return(prior)
}

.check_FUNL <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  if (distribution == "normal"){
    temp = prior[[1]]
    if (nrow(temp) != 7){
      stop("Normal Power model prior requires 7 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "FUNL Model [normal]"
    prior$parameters <- c("a","b","lm","ls","nm","ns","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) != 8){
      stop("Normal-NCV FUNL model prior requires 8 parameters.")
    }
  
    prior$model = "FUNL Model [normal-ncv]"
    prior$parameters <- c("a","b","lm","ls","nm","ns","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/FUNL specification is not presently available.")
  }
  prior$mean = .continuous_models[5] #give label in global varaible .continuous_models
  return(prior)
}


.check_exp3 <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
  
  if (distribution == "normal"){
    
    if (nrow(temp) != 4){
      stop("Normal Exponential-3  model prior requires 4 parameters.")
    }
    if (sum(temp[,4] > temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    
    prior$model = "Exponential-3 Model [normal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    if (nrow(temp) != 5){
      stop("Normal Exponential-3 model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[4,4] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Exponential-3 [normal-ncv]"
    prior$parameters <- c("a","b","c","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    temp = prior[[1]]
    if (nrow(temp) != 4){
      stop("Lognormal Exponential-3  model prior requires 4 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Exponential-3 Model [lognormal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  prior$mean = .continuous_models[2]
  temp <- prior$prior
  prior$prior = matrix(NA,nrow=nrow(temp)+1,5)
  prior$prior[1:2,] = temp[1:2,]
  prior$prior[3,]   = c(1,0,1,-100,100)
  prior$prior[4:nrow(prior$prior), ] = temp[3:nrow(temp),]
 
  cat("NOTE: Parameter 'c' added to prior list. It is not used in the analysis.\n")
  return(prior)
}

.check_polynomial <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
   if (sum(temp[,4]>temp[,5])> 0){
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  
  temp_p <- c("b0")
  for (ii in 2:(nrow(temp)-1)){
    temp_p <- c(temp_p,sprintf("b%s",ii-1))
  }
  
  if (distribution == "normal"){
    if (nrow(temp) < 3){
      stop("Normal Polynomial models require 3 or more parameters.")
    }
    prior$model = "Polynomial Model [normal]"
    prior$parameters <- c(temp_p,"log(sigma^2)")
    prior$degree = nrow(temp) - 2
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) < 4){
      stop("Normal-ncv polynomial models require 4 or more parameters.")
    }

    if (temp[nrow(temp)-1,4] < 0){
      stop("The prior on \\rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Polynomial Model [normal-ncv]"
    temp_p[length(temp_p)] = "rho"
    prior$parameters <- c(temp_p,"log(sigma^2)")
    prior$degree = nrow(temp) - 3
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/Polynomial specification is not presently available.")
  }
   
  prior$mean = .continuous_models[6]
  return(prior)
}

.check_FUNLhill <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
  if (sum(temp[,4]>temp[,5])> 0){
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  
  if (distribution == "normal"){
    
    if (nrow(temp) != 7){
      stop("Normal FUNL model prior requires 7 parameters.")
    }
  
    prior$model = "FUNL Model [normal]"
    prior$parameters <- c("a","b","LM","LD","NM","ND","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){

    if (nrow(temp) != 8){
      stop("Normal-NCV Hill model prior requires 8 parameters.")
    }
    
    if (temp[7,5] < 0){ #check rho
      stop("The prior on \rho (parameter 7) can not have a lower bound less than zero.")
    }
    prior$model = "FUNL Model [normal-ncv]"
    prior$parameters <- c("a","b","LM","LD","NM","ND","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/FUNL specification is not presently available.")
  }
  prior$mean = .continuous_models[5]
  return(prior)
}


#' @title create_dichotomous_prior  Given priorlist, a model, 
#'        and a distribution. Create a prior for a given analysis. 
#' @param prior First Prior
#' @param model Model to be used should be one of"hill","gamma","logistic","log-logistic","log-probit","multistage", "probit", "qlinear", or "weibull" 
#' @return new BMDprior list that can be used in a dichotomous fit. 
#' 
#' @examples 
#' plist<- create_prior_list(normprior(0,0.1,-100,100), # a
#'                           lnormprior(1,0.2,0,18))
#'  
#'  power_normal <- create_dichotomous_prior(plist,"logistic") 
#'
create_dichotomous_prior <- function(prior,model){
  
  if (!("BMDmodelprior" %in% class(prior))){
    stop("Prior is not of a 'BMDmodelprior' class. A probable solution is to 
          define the prior using function `create_prior_list`.")
  }
  if (!(model %in% .dichotomous_models  )){
    stop(cat("Model Type must be one of:", .dichotomous_models ,"\n"))
  }

  
  p = NA
  temp = prior[[1]]
  
  if (sum(temp[,4]>temp[,5])> 0){
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  
  if ("hill" == model){
    p = .check_d_hill(prior) 
  }
  
  if ("gamma" == model){
    p = .check_d_gamma(prior) 
  }
  
  if ("logistic" == model){
    p = .check_d_logistic(prior) 
  }
  
  if ("log-logistic" == model){
    p = .check_d_llogistic(prior) 
  }
  
  if ("multistage" == model){
    p = .check_d_multistage(prior) 
  }
  
  if ("probit" == model){
    p = .check_d_probit(prior) 
  }
  
  if ("log-probit" == model){
    p = .check_d_lprobit(prior) 
  }
  
  if ("qlinear" == model){
    p = .check_d_qlinear(prior) 
  }
  
  if ("weibull" == model){
    p = .check_d_weibull(prior) 
  }
  
  class(p)<- "BMD_Bayes_dichotomous_model"
  
  return(p)
}

.check_d_gamma   <- function(prior){
  
    temp = prior[[1]]
    
    if (nrow(temp) != 3){
      stop("Dichotomous Gamma  model prior requires 3 parameters.")
    }
    
    if (temp[2,4] < 0){
      stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
    }
    
    if (temp[3,4] < 0){
      stop("The prior on d (parameter 3) can not have a lower bound less than zero.")
    }
    
    prior$model = "Gamma Model [binomial]"
    prior$mean  = "gamma"
    prior$parameters <- c("logit(g)","a","b")
    return(prior)
}

.check_d_weibull <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 3){
    stop("Dichotomous Weibull  model prior requires 3 parameters.")
  }
  
  if (temp[2,4] < 0){
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }
  
  if (temp[3,4] < 0){
    stop("The prior on d (parameter 3) can not have a lower bound less than zero.")
  }
  

  prior$model = "Weibull Model [binomial]"
  prior$mean  = "weibull"
  prior$parameters <- c("logit(g)","a","b")
  return(prior)
}

.check_d_lprobit <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 3){
    stop("Dichotomous Log-Probit  model prior requires 3 parameters.")
  }

  if (temp[3,4] < 0){
    stop("The prior on b1 (parameter 3) can not have a lower bound less than zero.")
  }

  prior$model = "Log-Probit Model [binomial]"
  prior$mean  = "log-probit"
  prior$parameters <- c("logit(g)","b0","b1")
  return(prior)
}

.check_d_llogistic <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 3){
    stop("Dichotomous Log-Logistic model prior requires 3 parameters.")
  }
  

  if (temp[3,4] < 0){
    stop("The prior on b1 (parameter 3) can not have a lower bound less than zero.")
  }
  
  
  prior$model = "Log-Logistic Model [binomial]"
  prior$mean  = "log-logistic"
  prior$parameters <- c("logit(g)","b0","b1")
  return(prior)
}

.check_d_logistic <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 2){
    stop("Dichotomous logistic model prior requires 2 parameters.")
  }
  
  if (temp[2,4] < 0){
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }
  
  prior$model = "Logistic Model [binomial]"
  prior$mean  = "logistic"
  prior$parameters <- c("a","b")
  return(prior)
}

.check_d_probit     <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 2){
    stop("Dichotomous probit model prior requires 2 parameters.")
  }
  
  if (temp[2,4] < 0){
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }
  
  prior$model = "Probit Model [binomial]"
  prior$mean  = "probit"
  prior$parameters <- c("a","b")
  return(prior)
}

.check_d_qlinear    <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 2){
    stop("Dichotomous Quantal  model prior requires 2 parameters.")
  }
  
  if (temp[2,4] < 0){
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }

  prior$model = "Quantal Linear Model [binomial]"
  prior$mean  = "qlinear"
  prior$parameters <- c("logit(g)","b")
  return(prior)
}

.check_d_multistage <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) < 2 ){
    stop("Multistage model prior requires 2 or more parameters.")
  }
  
  names <- c("logit(g)")
  for (ii in 2:nrow(temp)){
    names <- c(names,sprintf("b%d",ii-1))
  }
  
  if (sum(temp[2:nrow(temp),4] < 0) > 0){
    stop("The prior on bx  can not have a lower bounds less than zero.")
  }
  rV <- list()
  rV$priors = temp
  rV$model = sprintf("Multistage-%d Model [binomial]",nrow(temp)-1)
  rV$mean  = "multistage"
  rV$degree = nrow(temp)-1
  rV$parameters <- names
  return(rV)
  
}

.check_d_hill <- function(prior){
  
  temp = prior[[1]]
  
  if (nrow(temp) != 4){
    stop("Dichotomous Hill  model prior requires 3 parameters.")
  }
  
  if (temp[4,4] < 0){
    stop("The prior on d (parameter 4) can not have a lower bound less than zero.")
  }
  prior$model = "Hill Model [binomial]"
  prior$mean  = "hill"
  prior$parameters <- c("logit(g)","b","c","d")
  return(prior)
}



