# Copyright 2020  NIEHS <matt.wheeler@nih.gov>
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.parse_prior <- function(prior) {
  rV <- list()
  rV$prior <- prior$prior

  temp_a <- regexpr("[[][a-zA-Z]+-*[a-zA-Z]+[]]", prior$model)
  start <- temp_a[1] + 1
  end <- start + attr(temp_a, "match.length") - 3
  if (temp_a == -1) {
    stop("Could not find a distribution for analysis.")
  }
  rV$distribution <- substr(prior$model, start, end)
  rV$model <- prior$mean

  class(rV) <- "BMD_Bayes_continuous_model"
  
  return(rV)
}

#' Build a list of default priors for model averaging
#'
#' @title ma_build_priors - Build a list of default priors for model averaging
#' @param Y Response matrix (optional) but preferred to fix variance prior
#' @param model_list Vector of models to use. See \code{\link{single_continuous_fit}} for details on the models. Overrides ma_type
#' @param distribution_list Vector of corresponding underlying distribution to use ("normal", "normal-ncv", "lognormal"). Defaults to normal-ncv
#' @param ma_type Models to average over ("original", "EFSA", or "all"). All corresponds to all implemented models; 
#' EFSA and original are described in \code{\link{ma_continuous_fit}} under the EFSA parameter.
#'
#' @return A list of configurations for the single models (priors, model type). Can be fed as the model_list to ma_continuous_fit
#'
#' @examples
#' #get EFSA model priors
#' prior_list <- ma_build_priors(ma_type = "EFSA")
#' 
ma_build_priors <- function(Y = NA, model_list = NA, distribution_list = NA, ma_type = "EFSA"){
  #build model/distribution lists by ma_type
  if(is.na(model_list[1])){
    if(ma_type == "EFSA"){
      model_list <- c(rep("exp-aerts", 2), rep("invexp-aerts",2), rep("hill-aerts", 2), rep("lognormal-aerts",2),
                      rep("gamma-efsa", 2), rep("LMS", 2), rep("probit-aerts", 2), rep("logistic-aerts",2))
      distribution_list <- rep(c("normal", "lognormal"), 8)
    }
    if(ma_type == "original"){
      model_list <- c(rep("hill", 2), rep("exp-3", 3), rep("exp-5", 3), rep("power", 2))
      distribution_list <- c(
        "normal", "normal-ncv", rep(c("normal", "normal-ncv", "lognormal"), 2),
        "normal", "normal-ncv"
      )
    }
    if(ma_type == "all"){
      model_list <- rep(.continuous_models, each=length(.continuous_distributions))
      distribution_list <- rep(.continuous_distributions, length(.continuous_models))
      badd_inds <- c(3, 12, 13:18)
      model_list <- model_list[-badd_inds]
      distribution_list <- distribution_list[-badd_inds]
    }
  }else{
    #check that the models are all valid
    if(!all(model_list %in% .continuous_models[-(5:6)])){
      stop('Please specify only the following model types: \n
            "hill","exp-3","exp-5","power", "exp-aerts", "invexp-aerts", "gamma-aerts", "invgamma-aerts", "hill-aerts",
            "lomax-aerts", "invlomax-aerts", "lognormal-aerts", "logskew-aerts", "invlogskew-aerts", "logistic-aerts", "probit-aerts", "LMS",
            "gamma-efsa"')
    }
    #if no distributions, default to normal-ncv
    if(is.na(distribution_list[1])){
      temp <- unique(model_list)
      if(length(temp) != length(model_list)){
        warning("Removing duplicate models. Please specify distribution_list to avoid this behavior.")
        model_list <- temp
      }
      distribution_list <- rep("normal-ncv", length(temp))
    }else{
      #check size compatibility
      if(length(model_list) != length(distribution_list)){
        stop("Please specify a distribution for each model.")
      }
      #check distribution validity
      if(!all(distribution_list %in% .continuous_distributions)){
        stop('Please specify only the following distribution types: "normal", "normal-ncv", "lognormal"')
      }
    }
  }
  #check for negative values if response provided
  if(!is.na(Y[1])){
    is_neg = .check_negative_response(Y)
    if(is_neg){
      tmpIdx <- which(distribution_list == "lognormal")
      model_list <- model_list[-tmpIdx]
      distribution_list <- distribution_list[-tmpIdx]
      if(identical(tmpIdx, integer(0))){
        warning("Negative response values were found in the data.  All lognormal
            models were removed from the analysis.")
      }
    }
  }
  prior_list <- list()
  for (ii in 1:length(model_list)) {
    PR <- .bayesian_prior_continuous_default(model_list[ii], distribution_list[ii], 2)
    # specify variance of last parameter to variance of response if response provided
    if(!is.na(Y[1])){
      if (distribution_list[ii] == "lognormal") {
        if (ncol(Y) > 1) {
          PR$priors[nrow(PR$priors), 2] <- log(mean(Y[, 3])^2)
        } else {
          PR$priors[nrow(PR$priors), 2] <- log(var(log(Y)))
        }
      } else {
        if (ncol(Y) > 1) {
          if (distribution_list[ii] == "normal") {
            PR$priors[nrow(PR$priors), 2] <- log(mean(Y[, 3])^2)
          } else {
            PR$priors[nrow(PR$priors), 2] <- log(abs(mean(Y[1, ])) / mean(Y[, 3])^2)
          }
        } else {
          if (distribution_list[ii] == "normal") {
            PR$priors[nrow(PR$priors), 2] <- log(var(Y))
          } else {
            PR$priors[nrow(PR$priors), 2] <- log(abs(mean(Y)) / var(Y))
          }
        }
      }
    }
    t_prior_result <- create_continuous_prior(PR, model_list[ii], distribution_list[ii], 2)
    # PR <- t_prior_result$prior
    prior_list[[ii]] <- t_prior_result#list(prior = PR, model_tye = model_list[ii], dist = distribution_list[ii])
  }
  return(prior_list)
}

#' @title create_continuous_prior  Given priorlist, a model,
#'        and a distribution. Create a prior for a given analysis.
#' @param prior_list First Prior
#' @param model Model to be used
#' one of \code{"hill", "exp-3", "exp-5", "power", "polynomial", "FUNL", "exp-aerts", "invexp-aerts", "gamma-aerts", "invgamma-aerts", "hill-aerts", 
#' "lomax-aerts", "invlomax-aerts", "lognormal-aerts", "logskew-aerts", "invlogskew-aerts", "logistic-aerts", "probit-aerts", "gamma-efsa"}
#' @param distribution - Normal "normal", Normal non-constant variance "normal-ncv", or
#'                       log-normal "lognormal"
#' @param deg - For polynomial models only, the degree of the polynomial.
#' @return new BMDprior list. This object is essentially a prior list constructed by
#' \code{create_prior_lists} with a model type and variance type.
#'
#' @examples
#' plist <- create_prior_list(
#'   normprior(0, 0.1, -100, 100), # a
#'   normprior(0, 1, -1e2, 1e2), # b
#'   lnormprior(1, 0.2, 0, 18), # k
#'   normprior(0, 1, -18, 18)
#' )
#'
#' power_normal <- create_continuous_prior(plist, "power", "normal")
#'
create_continuous_prior <- function(prior_list, model, distribution, deg = 2) {
  if (!("BMDmodelprior" %in% class(prior_list))) {
    stop("Prior is not of a 'BMDmodelprior' class. A probable solution is to
          define the prior using function `create_prior_list`.")
  }
  if (!(model %in% .continuous_models)) {
    stop(cat("Model Type must be one of:", .continuous_models, "\n"))
  }
  if (!(distribution %in% .continuous_distributions)) {
    stop(cat("Distribution must be one of:", .continuous_distributions, "\n"))
  }
  temp <- floor(deg)

  if (deg < 1) {
    stop("Polynomial degree must be greater than or equal to 1.")
  }
  if (temp != deg) {
    stop("Polynomial degree must be an integer.")
  }

  p <- NA
  if(grepl("aerts", model)){
    if(model %in% c("exp-aerts", "invexp-aerts", "hill-aerts", "lognormal-aerts")){
      p <- .check_4param(prior_list, distribution)
    }
    if(model %in% c("logistic-aerts", "probit-aerts")){
      p <- .check_4param_sigmoidal(prior_list, distribution)
    }
    if(model %in% c("gamma-aerts", "invgamma-aerts", "lomax-aerts", "invlomax-aerts",
                    "logskew-aerts", "invlogskew-aerts")){
      p <- .check_5param(prior_list, distribution)
    }
    if(model %in% c("gamma-aerts","invgamma-aerts")){
      p<- .check_5param_gamma(prior_list,distribution)
    }
    p$mean <- model
  }else{
    if("LMS" == model){
      p <- .check_4param(prior_list, distribution)
      p$mean <- model
    }
    if("gamma-efsa" == model){
      p <- .check_4param_gamma(prior_list, distribution)
      p$mean <- model
    }
    if ("hill" == model) {
      p <- .check_hill(prior_list, distribution)
    }
    
    if ("FUNL" == model) {
      p <- .check_FUNL(prior_list, distribution)
    }
    
    if ("exp-5" == model) {
      p <- .check_exp5(prior_list, distribution)
    }
    
    if ("exp-3" == model) {
      p <- .check_exp3(prior_list, distribution)
    }
    
    if ("polynomial" == model) {
      p <- .check_polynomial(prior_list, distribution)
    }
    
    if ("power" == model) {
      p <- .check_power(prior_list, distribution)
    }
  }

  class(p) <- "BMD_Bayes_continuous_model"

  return(p)
}

.check_hill <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  if (distribution == "normal") {
    temp <- prior[[1]]
    if (nrow(temp) != 5) {
      stop("Normal Hill model prior requires 5 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model <- "Hill Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }

  if (distribution == "normal-ncv") {
    temp <- prior[[1]]
    if (nrow(temp) != 6) {
      stop("Normal-NCV Hill model prior requires 6 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5, 4] < 0) {
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Hill Model [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "rho", "log(sigma^2)")
  }

  if (distribution == "lognormal") {
    stop("Log-Normal/Hill specification is not presently available.")
  }
  prior$mean <- .continuous_models[1]
  return(prior)
}

.check_4param_gamma <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  if (temp[2,4] <= 0){
    stop("The lower bound on b must be positive (for numerical stability).")
  }
  if(temp[4,4] < 0.2){
    stop("The lower bound on d must be greater than or equal to 0.2 for numerical reasons. ")
  }
  
  if (distribution == "normal") {
    if (nrow(temp) != 5) {
      stop("Normal Aerts model prior requires 5 parameters.")
    }
    prior$model <- "Aerts Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }
  
  if (distribution == "normal-ncv") {
    if (nrow(temp) != 6) {
      stop("Normal-NCV Aerts model prior requires 6 parameters.")
    }
    if (temp[5, 4] < 0) {
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Aerts Model [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "rho", "log(sigma^2)")
  }
  
  if (distribution == "lognormal") {
    if (nrow(temp) != 5) {
      stop("Lognormal Aerts model prior requires 5 parameters.")
    }
    prior$model <- "Aerts Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")  }
  return(prior)
}


.check_4param <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  if (temp[2,4] <= 0){
    stop("The lower bound on b must be positive (for numerical stability).")
  }
  if(temp[4,4] < 0){
    stop("The lower bound on d must be non-negative.")
  }
  if (distribution == "normal") {
    if (nrow(temp) != 5) {
      stop("Normal Aerts model prior requires 5 parameters.")
    }
    prior$model <- "Aerts Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }
  
  if (distribution == "normal-ncv") {
    if (nrow(temp) != 6) {
      stop("Normal-NCV Aerts model prior requires 6 parameters.")
    }
    if (temp[5, 4] < 0) {
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Aerts Model [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "rho", "log(sigma^2)")
  }
  
  if (distribution == "lognormal") {
    if (nrow(temp) != 5) {
      stop("Lognormal Aerts model prior requires 5 parameters.")
    }
    prior$model <- "Aerts Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")  }
  return(prior)
}

.check_4param_sigmoidal <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  if(temp[4,4] < 0){
    stop("The lower bound on d (the exponent) must be non-negative.")
  }
  if (distribution == "normal") {
    if (nrow(temp) != 5) {
      stop("Normal Aerts model prior requires 5 parameters.")
    }
    prior$model <- "Aerts Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }
  
  if (distribution == "normal-ncv") {
    if (nrow(temp) != 6) {
      stop("Normal-NCV Aerts model prior requires 6 parameters.")
    }
    if (temp[5, 4] < 0) {
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Aerts Model [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "rho", "log(sigma^2)")
  }
  
  if (distribution == "lognormal") {
    if (nrow(temp) != 5) {
      stop("Lognormal Aerts model prior requires 5 parameters.")
    }
    prior$model <- "Aerts Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")  }
  return(prior)
}

.check_5param <- function(prior, distribution) {
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  if (temp[2,4] <= 0  | temp[5,4] <= 0){
    stop("The lower bounds on b and xi must be positive (for numerical stability).")
  }
  if(temp[4,4] < 0){
    stop("The lower bound on d must be non-negative.")
  }
  if(temp[2,4] < 0.1){
    warning("Lower bound on b < 0.1 may cause numerical instability")
  }
  # check if the normal distribution is correctly specified
  if (distribution == "normal") {
    if (nrow(temp) != 6) {
      stop("Normal Aerts model prior requires 6 parameters.")
    }
    prior$model <- "Aerts Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "xi", "log(sigma^2)")
  }
  
  if (distribution == "normal-ncv") {
    temp <- prior[[1]]
    if (nrow(temp) != 7) {
      stop("Normal-NCV Aerts model prior requires 7 parameters.")
    }
    if (temp[6, 4] < 0) {
      stop("The prior on \rho (parameter 6) can not have a lower bound less than zero.")
    }
    prior$model <- "Aerts Model [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "xi", "rho", "log(sigma^2)")
  }
  
  if (distribution == "lognormal") {
    temp <- prior[[1]]
    if (nrow(temp) != 6) {
      stop("Lognormal Aerts model prior requires 6 parameters.")
    }
    prior$model <- "Aerts Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "xi", "log(sigma^2)")
  }
  prior$mean <- .continuous_models[1]
  return(prior)
}


.check_5param_gamma <- function(prior, distribution) {
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  if (temp[2,4] <= 0  | temp[5,4] <= 0){
    stop("The lower bounds on b and xi must be positive (for numerical stability).")
  }
  if(temp[4,4] < 0){
    stop("The lower bound on d must be non-negative.")
  }
   if(temp[5,4] < 0.2){
    stop("The lower bound on xi must be greater than or equal to 0.2.")
  }
  if(temp[2,4] < 0.1){
    warning("Lower bound on b < 0.1 may cause numerical instability")
  }
  # check if the normal distribution is correctly specified
  if (distribution == "normal") {
    if (nrow(temp) != 6) {
      stop("Normal Aerts model prior requires 6 parameters.")
    }
    prior$model <- "Aerts Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "xi", "log(sigma^2)")
  }
  
  if (distribution == "normal-ncv") {
    temp <- prior[[1]]
    if (nrow(temp) != 7) {
      stop("Normal-NCV Aerts model prior requires 7 parameters.")
    }
    if (temp[6, 4] < 0) {
      stop("The prior on \rho (parameter 6) can not have a lower bound less than zero.")
    }
    prior$model <- "Aerts Model [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "xi", "rho", "log(sigma^2)")
  }
  
  if (distribution == "lognormal") {
    temp <- prior[[1]]
    if (nrow(temp) != 6) {
      stop("Lognormal Aerts model prior requires 6 parameters.")
    }
    prior$model <- "Aerts Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "xi", "log(sigma^2)")
  }
  prior$mean <- .continuous_models[1]
  return(prior)
}

.check_exp5 <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]
  if (distribution == "normal") {
    if (nrow(temp) != 5) {
      stop("Normal Exponential-5  model prior requires 5 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }

    prior$model <- "Exponential-5 Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }

  if (distribution == "normal-ncv") {
    if (nrow(temp) != 6) {
      stop("Normal Exponential-5 model prior requires 6 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5, 4] < 0) {
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Exponential-5 [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "rho", "log(sigma^2)")
  }

  if (distribution == "lognormal") {
    temp <- prior[[1]]
    if (nrow(temp) != 5) {
      stop("Lognormal Exponential-5  model prior requires 5 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model <- "Exponential-5 Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }
  prior$mean <- .continuous_models[3]
  return(prior)
}

.check_power <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  if (distribution == "normal") {
    temp <- prior[[1]]
    if (nrow(temp) != 4) {
      stop("Normal Power model prior requires 5 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[3, 4] < 0) {
      stop("The power parameter d (parameter 3) can not have a lower bound less than zero.")
    }
    prior$model <- "Power Model [normal]"
    prior$parameters <- c("a", "b", "d", "log(sigma^2)")
  }

  if (distribution == "normal-ncv") {
    temp <- prior[[1]]
    if (nrow(temp) != 5) {
      stop("Normal-NCV Power model prior requires 5 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[3, 4] < 0) {
      stop("The power parameter d (parameter 3) can not have a lower bound less than zero.")
    }

    if (temp[4, 4] < 0) {
      stop("The prior on \rho (parameter 4) can not have a lower bound less than zero.")
    }

    prior$model <- "Power Model [normal-ncv]"
    prior$parameters <- c("a", "b", "d", "rho", "log(sigma^2)")
  }

  if (distribution == "lognormal") {
    stop("Log-Normal/Power specification is not presently available.")
  }
  prior$mean <- .continuous_models[4]
  return(prior)
}

.check_FUNL <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  if (distribution == "normal") {
    temp <- prior[[1]]
    if (nrow(temp) != 7) {
      stop("Normal Power model prior requires 7 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model <- "FUNL Model [normal]"
    prior$parameters <- c("a", "b", "lm", "ls", "nm", "ns", "log(sigma^2)")
  }

  if (distribution == "normal-ncv") {
    temp <- prior[[1]]
    if (nrow(temp) != 8) {
      stop("Normal-NCV FUNL model prior requires 8 parameters.")
    }

    prior$model <- "FUNL Model [normal-ncv]"
    prior$parameters <- c("a", "b", "lm", "ls", "nm", "ns", "rho", "log(sigma^2)")
  }

  if (distribution == "lognormal") {
    stop("Log-Normal/FUNL specification is not presently available.")
  }
  prior$mean <- .continuous_models[5] # give label in global varaible .continuous_models
  return(prior)
}


.check_exp3 <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]

  if (distribution == "normal") {
    if (nrow(temp) != 4) {
      stop("Normal Exponential-3  model prior requires 4 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }

    prior$model <- "Exponential-3 Model [normal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }

  if (distribution == "normal-ncv") {
    if (nrow(temp) != 5) {
      stop("Normal Exponential-3 model prior requires 5 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[4, 4] < 0) {
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Exponential-3 [normal-ncv]"
    prior$parameters <- c("a", "b", "c", "d", "rho", "log(sigma^2)")
  }

  if (distribution == "lognormal") {
    temp <- prior[[1]]
    if (nrow(temp) != 4) {
      stop("Lognormal Exponential-3  model prior requires 4 parameters.")
    }
    if (sum(temp[, 4] > temp[, 5]) > 0) {
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model <- "Exponential-3 Model [lognormal]"
    prior$parameters <- c("a", "b", "c", "d", "log(sigma^2)")
  }

  prior$mean <- .continuous_models[2]
  temp <- prior$prior
  prior$prior <- matrix(NA, nrow = nrow(temp) + 1, 5)
  prior$prior[1:2, ] <- temp[1:2, ]
  prior$prior[3, ] <- c(1, 0, 1, -100, 100)
  prior$prior[4:nrow(prior$prior), ] <- temp[3:nrow(temp), ]

  cat("NOTE: Parameter 'c' added to prior list. It is not used in the analysis.\n")
  return(prior)
}

.check_polynomial <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }

  temp_p <- c("b0")
  for (ii in 2:(nrow(temp) - 1)) {
    temp_p <- c(temp_p, sprintf("b%s", ii - 1))
  }

  if (distribution == "normal") {
    if (nrow(temp) < 3) {
      stop("Normal Polynomial models require 3 or more parameters.")
    }
    prior$model <- "Polynomial Model [normal]"
    prior$parameters <- c(temp_p, "log(sigma^2)")
    prior$degree <- nrow(temp) - 2
  }

  if (distribution == "normal-ncv") {
    temp <- prior[[1]]
    if (nrow(temp) < 4) {
      stop("Normal-ncv polynomial models require 4 or more parameters.")
    }

    if (temp[nrow(temp) - 1, 4] < 0) {
      stop("The prior on \\rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model <- "Polynomial Model [normal-ncv]"
    temp_p[length(temp_p)] <- "rho"
    prior$parameters <- c(temp_p, "log(sigma^2)")
    prior$degree <- nrow(temp) - 3
  }

  if (distribution == "lognormal") {
    stop("Log-Normal/Polynomial specification is not presently available.")
  }

  prior$mean <- .continuous_models[6]
  return(prior)
}

.check_FUNLhill <- function(prior, distribution) {
  # check if the normal distribution is correctly specified
  temp <- prior[[1]]
  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }

  if (distribution == "normal") {
    if (nrow(temp) != 7) {
      stop("Normal FUNL model prior requires 7 parameters.")
    }

    prior$model <- "FUNL Model [normal]"
    prior$parameters <- c("a", "b", "LM", "LD", "NM", "ND", "log(sigma^2)")
  }

  if (distribution == "normal-ncv") {
    if (nrow(temp) != 8) {
      stop("Normal-NCV Hill model prior requires 8 parameters.")
    }

    if (temp[7, 5] < 0) { # check rho
      stop("The prior on \rho (parameter 7) can not have a lower bound less than zero.")
    }
    prior$model <- "FUNL Model [normal-ncv]"
    prior$parameters <- c("a", "b", "LM", "LD", "NM", "ND", "rho", "log(sigma^2)")
  }

  if (distribution == "lognormal") {
    stop("Log-Normal/FUNL specification is not presently available.")
  }
  prior$mean <- .continuous_models[5]
  return(prior)
}


#' @title create_dichotomous_prior  Given priorlist, a model,
#'        and a distribution. Create a prior for a given analysis.
#' @param prior First Prior
#' @param model Model to be used should be one of"hill","gamma","logistic","log-logistic","log-probit","multistage", "probit", "qlinear", or "weibull"
#' @return new BMDprior list that can be used in a dichotomous fit.
#'
#' @examples
#' plist <- create_prior_list(
#'   normprior(0, 0.1, -100, 100), # a
#'   lnormprior(1, 0.2, 0, 18)
#' )
#'
#' power_normal <- create_dichotomous_prior(plist, "logistic")
#'
create_dichotomous_prior <- function(prior, model) {
  if (!("BMDmodelprior" %in% class(prior))) {
    stop("Prior is not of a 'BMDmodelprior' class. A probable solution is to
          define the prior using function `create_prior_list`.")
  }
  if (!(model %in% .dichotomous_models)) {
    stop(cat("Model Type must be one of:", .dichotomous_models, "\n"))
  }


  p <- NA
  temp <- prior[[1]]

  if (sum(temp[, 4] > temp[, 5]) > 0) {
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }

  if ("hill" == model) {
    p <- .check_d_hill(prior)
  }

  if ("gamma" == model) {
    p <- .check_d_gamma(prior)
  }

  if ("logistic" == model) {
    p <- .check_d_logistic(prior)
  }

  if ("log-logistic" == model) {
    p <- .check_d_llogistic(prior)
  }

  if ("multistage" == model) {
    p <- .check_d_multistage(prior)
  }

  if ("probit" == model) {
    p <- .check_d_probit(prior)
  }

  if ("log-probit" == model) {
    p <- .check_d_lprobit(prior)
  }

  if ("qlinear" == model) {
    p <- .check_d_qlinear(prior)
  }

  if ("weibull" == model) {
    p <- .check_d_weibull(prior)
  }

  class(p) <- "BMD_Bayes_dichotomous_model"

  return(p)
}

.check_d_gamma <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 3) {
    stop("Dichotomous Gamma  model prior requires 3 parameters.")
  }

  if (temp[2, 4] < 0) {
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }

  if (temp[3, 4] < 0) {
    stop("The prior on d (parameter 3) can not have a lower bound less than zero.")
  }

  prior$model <- "Gamma Model [binomial]"
  prior$mean <- "gamma"
  prior$parameters <- c("logit(g)", "a", "b")
  return(prior)
}

.check_d_weibull <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 3) {
    stop("Dichotomous Weibull  model prior requires 3 parameters.")
  }

  if (temp[2, 4] < 0) {
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }

  if (temp[3, 4] < 0) {
    stop("The prior on d (parameter 3) can not have a lower bound less than zero.")
  }


  prior$model <- "Weibull Model [binomial]"
  prior$mean <- "weibull"
  prior$parameters <- c("logit(g)", "a", "b")
  return(prior)
}

.check_d_lprobit <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 3) {
    stop("Dichotomous Log-Probit  model prior requires 3 parameters.")
  }

  if (temp[3, 4] < 0) {
    stop("The prior on b1 (parameter 3) can not have a lower bound less than zero.")
  }

  prior$model <- "Log-Probit Model [binomial]"
  prior$mean <- "log-probit"
  prior$parameters <- c("logit(g)", "b0", "b1")
  return(prior)
}

.check_d_llogistic <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 3) {
    stop("Dichotomous Log-Logistic model prior requires 3 parameters.")
  }


  if (temp[3, 4] < 0) {
    stop("The prior on b1 (parameter 3) can not have a lower bound less than zero.")
  }


  prior$model <- "Log-Logistic Model [binomial]"
  prior$mean <- "log-logistic"
  prior$parameters <- c("logit(g)", "b0", "b1")
  return(prior)
}

.check_d_logistic <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 2) {
    stop("Dichotomous logistic model prior requires 2 parameters.")
  }

  if (temp[2, 4] < 0) {
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }

  prior$model <- "Logistic Model [binomial]"
  prior$mean <- "logistic"
  prior$parameters <- c("a", "b")
  return(prior)
}

.check_d_probit <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 2) {
    stop("Dichotomous probit model prior requires 2 parameters.")
  }

  if (temp[2, 4] < 0) {
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }

  prior$model <- "Probit Model [binomial]"
  prior$mean <- "probit"
  prior$parameters <- c("a", "b")
  return(prior)
}

.check_d_qlinear <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 2) {
    stop("Dichotomous Quantal  model prior requires 2 parameters.")
  }

  if (temp[2, 4] < 0) {
    stop("The prior on b (parameter 2) can not have a lower bound less than zero.")
  }

  prior$model <- "Quantal Linear Model [binomial]"
  prior$mean <- "qlinear"
  prior$parameters <- c("logit(g)", "b")
  return(prior)
}

.check_d_multistage <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) < 2) {
    stop("Multistage model prior requires 2 or more parameters.")
  }

  names <- c("logit(g)")
  for (ii in 2:nrow(temp)) {
    names <- c(names, sprintf("b%d", ii - 1))
  }

  if (sum(temp[2:nrow(temp), 4] < 0) > 0) {
    stop("The prior on bx  can not have a lower bounds less than zero.")
  }
  rV <- list()
  rV$priors <- temp
  rV$model <- sprintf("Multistage-%d Model [binomial]", nrow(temp) - 1)
  rV$mean <- "multistage"
  rV$degree <- nrow(temp) - 1
  rV$parameters <- names
  return(rV)
}

.check_d_hill <- function(prior) {
  temp <- prior[[1]]

  if (nrow(temp) != 4) {
    stop("Dichotomous Hill  model prior requires 3 parameters.")
  }

  if (temp[4, 4] < 0) {
    stop("The prior on d (parameter 4) can not have a lower bound less than zero.")
  }
  prior$model <- "Hill Model [binomial]"
  prior$mean <- "hill"
  prior$parameters <- c("logit(g)", "b", "c", "d")
  return(prior)
}
