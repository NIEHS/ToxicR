#################################################
# Prior File
#
#################################################
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


#' Specify a normal prior for a ToxicR Bayesian model fit.
#' @title normprior - create a normal prior object
#' @param mean mean of the prior distribution.
#' @param sd sd of the prior distribution.
#' @param lb lower bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully normal prior, make lb small relative to the mean/sd.
#' @param ub  Upper bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully normal prior, make ub large relative to the mean/sd.
#' @return a normal prior model object.  This object is essentially a vector with
#' the first element as 1 (for normal), the second element the mean, the third
#' element the standard deviation, the fourth and fifth elements the lower and upper bounds, respectively.
#'
#' @examples
#' # Normal Prior with mean 0,sd-1
#' normprior(mean = 0, sd = 1, lb = -1e4, ub = 1e4)
#'
#' # Truncated Normal prior, Truncated below at 0
#' normprior(mean = 0, sd = 1, lb = 0, ub = 1e4)
#'
normprior <- function(mean = 0, sd = 1, lb = -100, ub = 100) {
  if (ub < lb) {
    stop("Upper Bound must be greater than lower bound")
  }
  retValue <- matrix(c(1, mean, sd, lb, ub), ncol = 5)
  class(retValue) <- "BMDprior"
  return(retValue)
}

#' Specify a log-normal prior for a ToxicR Bayesian model fit.
#' @title lnormprior - create a lognormal prior.
#' @param mean log-mean of the prior distribution.
#' @param sd log-sd of the prior distribution.
#' @param lb lower bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully normal prior, make lb small relative to the mean/sd.
#' @param ub  Upper bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully normal prior, make ub large relative to the mean/sd.
#' @return a normal prior model object
#' This object is essentially a vector with
#' the first element as 2 (for log-normal), the second element the mean, the third
#' element the log-standard deviation, the fourth and fifth elements the lower and upper bounds, respectively.
#' @examples
#' # Log-Normal Prior with mean 0,sd-1
#' lnormprior(mean = 0, sd = 1, lb = -1e4, ub = 1e4)
#'
#' # Truncated Log-Normal prior, Truncated below at 1
#' lnormprior(mean = 0, sd = 1, lb = 1, ub = 1e4)
#'
lnormprior <- function(mean = 0, sd = 1, lb = -100, ub = 100) {
  if (lb < 0) {
    lb <- 0
  }

  if (ub < lb) {
    stop("Upper Bound must be greater than lower bound")
  }

  retValue <- matrix(c(2, mean, sd, lb, ub), ncol = 5)
  class(retValue) <- "BMDprior"
  return(retValue)
}

#' Specify a Cauchy prior for a ToxicR Bayesian model fit
#' @title cauchyprior - create a Cauchy prior.
#' @param mean mean of the prior distribution.
#' @param sd sd of the prior distribution.
#' @param lb lower bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully Cauchy prior, make lb small relative to the mean/sd.
#' @param ub  Upper bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully Cauchy prior, make ub large relative to the mean/sd.
#' @return a normal prior model object.  This object is essentially a vector with
#' the first element as 4 (for Cauchy), the second element the mean, the third
#' element the standard deviation, the fourth and fifth elements the lower and upper bounds, respectively.
#' @examples
#' # Cauchy Prior with mean 0,sd-1
#' cauchyprior(mean = 0, sd = 1, lb = -1e4, ub = 1e4)
#'
#' # Half Cauchy prior, Truncated below at 0
#' cauchyprior(mean = 0, sd = 1, lb = 0, ub = 1e4)
cauchyprior <- function(mean = 0, sd = 1, lb = -100, ub = 100) {
  if (ub < lb) {
    stop("Upper Bound must be greater than lower bound")
  }
  retValue <- matrix(c(3, mean, sd, lb, ub), ncol = 5)
  class(retValue) <- "BMDprior"
  return(retValue)
}

#' Specify a Gamma prior for a ToxicR Bayesian model fit
#' @title gammaprior - create a Gamma prior.
#' @param mean mean of the prior distribution.
#' @param sd sd of the prior distribution.
#' @param lb lower bound on the distribution. Necessary for the optimization algorithms, 
#' must be >= 0
#' @param ub  Upper bound on the distribution. Necessary for the optimization algorithms,
#' To make sure it is a fully Gamma prior, make ub large relative to the mean/sd.
#' @return a normal prior model object.  This object is essentially a vector with
#' the first element as 4 (for Gamma), the second element the mean, the third
#' element the standard deviation, the fourth and fifth elements the lower and upper bounds, respectively.
#' @examples
#' # Gamma Prior with mean 1,sd 1, and lower bound 0.1
#' gammaprior(mean = 1, sd = 1, lb = 0.1, ub = 1e4)
gammaprior <- function(mean = 1, sd = 1, lb = 0, ub = 100){
  if (ub < lb) {
    stop("Upper Bound must be greater than lower bound")
  }
  if(lb < 0){
    stop("Gamma distribution must be positive! Use 0 or greater")
  }
  retValue <- matrix(c(4, mean, sd, lb, ub), ncol = 5)
  class(retValue) <- "BMDprior"
  return(retValue)
}

#' Specify a modified-PERT (essentially a continuous version of a triangular distribution) prior for a ToxicR Bayesian model fit
#' @title pertprior - create a PERT prior.
#' @param mode mode of the prior distribution.
#' @param shape shape of the modified PERT prior distribution (often denoted as gamma), optional
#' @param lb lower bound on the distribution. Necessary for the optimization algorithms.
#' @param ub  Upper bound on the distribution. Necessary for the optimization algorithms,
#' @return a normal prior model object.  This object is essentially a vector with
#' the first element as 5 (for PERT), the second element the mean, the third
#' element the standard deviation, the fourth and fifth elements the lower and upper bounds, respectively.
#' @examples
#' # standard PERT Prior with mode 3, lower bound 1, upper bound 5
#' pertprior(mode = 3, lb = 1, ub = 5)
#' #modified PERT Prior with smaller shape parameter
#' pertprior(mode = 3, lb = 1, ub = 5, shape = 3)
pertprior <- function(mode = 1, shape = 4, lb = 0.1, ub = 100){
  if (ub < lb) {
    stop("Upper Bound must be greater than lower bound")
  }
  if(shape < 0){
    stop("Shape should be positive")
  }
  mean <- (lb + shape * mode + ub) / (shape + 2)
  sd <- sqrt((mean - lb) * (ub - mean) / (shape + 3))
  retValue <- matrix(c(5, mean, sd, lb, ub), ncol = 5)
  class(retValue) <- "BMDprior"
  return(retValue)
}





#' @title create_prior_lists .. Given priors
#'        created using the ToxicR prior functions, create a list of priors
#'        for a model.
#' @param x1 First Prior
#' @param x2 Second Prior
#' @param ... Aditional arguments
#' @return new BMDprior list. This object is essentailly a matrix where each
#' row is  an element defined by a prior object (e.g., normprior or lnormprior).
#'
#' @examples
#' plist <- create_prior_list(
#'   normprior(0, 0.1, -100, 100), # a
#'   normprior(0, 1, -1e2, 1e2), # b
#'   lnormprior(1, 0.2, 0, 18), # k
#'   normprior(0, 1, -18, 18)
#' )
#'
create_prior_list <- function(x1, x2, ...) {
  cl <- match.call()
  mf <- as.list(match.call(expand.dots = TRUE))[-1]

  X <- matrix(0, nrow = length(mf), ncol = 5)
  for (ii in 1:length(mf)) {
    X[ii, ] <- eval(mf[[ii]])
  }
  Y <- list()
  Y[[1]] <- X
  names(Y) <- c("priors")
  class(Y) <- "BMDmodelprior"
  return(Y)
}


.combine_prior_lists <- function(p1, p2) {
  if (as.character(class(p1)) == "BMDprior") {
    x1 <- as.matrix(p1[, , drop = F])
  } else {
    x1 <- p1[[1]]
  }

  if (as.character(class(p2)) == "BMDprior") {
    x2 <- as.matrix(p2[, , drop = F])
  } else {
    x2 <- p2[[1]]
  }

  retval <- list(priors = rbind(x1, x2))

  class(retval) <- "BMDmodelprior"
  return(retval)
}

.print.BMD_Bayes_model <- function(x, ...) {
  priors <- x
  X <- priors[[1]]
  if (!is.null(priors$model)) {
    cat(priors$model, " Parameter Priors\n")
  } else {
    cat("Model Parameter Priors\n ")
  }

  cat("------------------------------------------------------------------------\n")
  for (ii in 1:nrow(X)) {
    V <- X[ii, ]
    if (!is.null(priors$parameters)) {
      temp_text <- sprintf("Prior [%s]:", priors$parameters[ii])
    } else {
      temp_text <- "Prior: "
    }
    if (V[1] == 1) {
      cat(sprintf(
        "%sNormal(mu = %1.2f, sd = %1.3f) 1[%1.2f,%1.2f]\n", temp_text, V[2],
        V[3], V[4], V[5]
      ))
    }
    if (V[1] == 2) {
      cat(sprintf(
        "%sLog-Normal(log-mu = %1.2f, log-sd = %1.3f) 1[%1.2f,%1.2f]\n", temp_text, V[2],
        V[3], V[4], V[5]
      ))
    }
    if (V[1] == 3) {
      cat(sprintf(
        "%sCauchy(mu = %1.2f, sd = %1.3f) 1[%1.2f,%1.2f]\n", temp_text, V[2],
        V[3], V[4], V[5]
      ))
    }
    if (V[1] == 4) {
      cat(sprintf(
        "%sGamma(mu = %1.2f, sd = %1.3f) 1[%1.2f,%1.2f]\n", temp_text, V[2],
        V[3], V[4], V[5]
      ))
    }
    if (V[1] == 5) {
      cat(sprintf(
        "%sPERT(mu = %1.2f, sd = %1.3f) 1[%1.2f,%1.2f]\n", temp_text, V[2],
        V[3], V[4], V[5]
      ))
    }
  }
}



.bayesian_prior_continuous_default <- function(model, distribution, degree = 2) {
  variance <- distribution
  # dmodel <- which(model == c("hill", "exp-3", "exp-5", "power", "polynomial"))
  dmodel <- which(model == .continuous_models[-5]) #remove FUNL
  dvariance <- which(variance == c("normal", "normal-ncv", "lognormal"))

  if (dmodel == 1) { # Hill
    if (dvariance == 1) {
      prior <- create_prior_list(
        normprior(1, 1, -100, 100),
        normprior(0, 1, -100, 100),
        lnormprior(0, 2, 0, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 1, -30, 30)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        normprior(1, 1, -100, 100),
        normprior(0, 1, -100, 100),
        lnormprior(0, 2, 0, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0, 1, 0, 100),
        normprior(-3, 1, -18, 18)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        normprior(1, 1, -100, 100),
        normprior(0, 2, -100, 100), # normprior(1,2,-18,18),
        lnormprior(0, 1, 0, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 1, -18, 18)
      )
    }
  } else if (dmodel == 2) { # Exponential-3
    if (dvariance == 1) {
      prior <- create_prior_list(
        normprior(0, 1, -100, 100),
        lnormprior(0, 2, 0, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 2, -18, 18)
      )
    } else if (dvariance == 2) { # NonConstant Normal Prior
      prior <- create_prior_list(
        normprior(0, 1, -100, 100),
        lnormprior(0, 2, 0, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0, 0.5, 0, 18),
        normprior(0, 1, -30, 30)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        lnormprior(0, 1, 0, 100),
        lnormprior(0, 2, 0, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 1, -18, 18)
      )
    }
  } else if (dmodel == 3) { # Exp-5
    if (dvariance == 1 || dvariance == 3) {
      prior <- create_prior_list(
        lnormprior(0, 1, 0, 100),
        normprior(0, 1, -100, 100),
        normprior(0, 1, -100, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        lnormprior(0, 1, 0, 100),
        normprior(0, 1, -100, 100),
        normprior(0, 1, -100, 100),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0, 0.5, 0, 18),
        normprior(0, 1, -30, 30)
      )
    }
  } else if (dmodel == 4) { # Power
    if (dvariance == 1) {
      prior <- create_prior_list(
        normprior(0, 1, -100, 100),
        normprior(0, 1, -1e2, 1e2),
        lnormprior(log(1.6), 0.4214036, 0, 40),
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        normprior(0, 1, -100, 100), # a
        normprior(0, 1, -1e2, 1e2), # b
        lnormprior(log(1.6), 0.4214036, 0, 40), # k
        lnormprior(0, 0.5, 0, 18),
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 3) {
      stop("Power-Log-normal models are not allowed. Please choose normal or normal non-constant variance. \n")
    }
  } else if (dmodel == 5) { # Polynomial
    cat("WARNING: Polynomial models may provide unstable estimates because of possible non-monotone behavior.\n")
    if (dvariance == 1) {
      prior <- create_prior_list(normprior(0, 5, -100, 100))
      for (ii in 1:degree) {
        prior <- .combine_prior_lists(
          prior,
          normprior(0, 5, -100, 100)
        )
      }
      prior <- .combine_prior_lists(prior, create_prior_list(normprior(0, 1, -18, 18)))
    } else if (dvariance == 2) {
      prior <- create_prior_list(normprior(0, 5, -100, 100))
      for (ii in 1:degree) {
        prior <- .combine_prior_lists(
          prior,
          normprior(0, 5, -100, 100)
        )
      }
      prior <- .combine_prior_lists(
        prior,
        create_prior_list(
          lnormprior(0, 1, 0, 100),
          normprior(0, 1, -18, 18)
        )
      )
    } else if (dvariance == 3) {
      stop("Polynomial-Log-normal models are not allowed. Please choose normal or normal non-constant variance. \n")
    }
  } else if (dmodel %in% c(6:7, 10, 13,18)){ #4 parameter Aerts models
    if (dvariance == 1) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(0, 3, 1e-6, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 3, -30, 30)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(0, 3, 1e-6, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0, 2, 0, 100),
        normprior(-2, 3, -30, 30)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(0, 2, 1e-6, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 3, -30, 30)
      )
    }
  } else if (dmodel %in% c(16,17)){ #logistic/probit
    if (dvariance == 1) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        normprior(0, 3, -100, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 3, -30, 30)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        normprior(0, 3, -100, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0, 2, 0, 100),
        normprior(-2, 3, -30, 30)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        normprior(0, 3, -100, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        normprior(0, 3, -30, 30)
      )
    }
  } else if (dmodel %in% c(11:12,14:15)){ #5 parameter Aerts models
    if (dvariance == 1) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 5, 0.1, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0.2, 0.5, 1e-6, 18),
        normprior(0, 3, -30, 30)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 5, 0.1, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0.2, 0.5, 1e-6, 18),
        lnormprior(0, 2, 0, 100),
        normprior(-2, 3, -30, 30)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 5, 0.1, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0.2, 0.5, 1e-6, 18),
        normprior(0, 3, -30, 30)
      )
    }
    if(dmodel == 8){
      prior$priors[2,] <- lnormprior(1,2,0.2,20)
    }
  } else if (dmodel %in% c(8:9)){ #5 parameter Aerts models Gamma
    if (dvariance == 1) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 5, 0.1, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0.2, 0.5, 0.2, 18),
        normprior(0, 3, -30, 30)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 5, 0.1, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0.2, 0.5, 0.2, 18),
        lnormprior(0, 2, 0, 100),
        normprior(-2, 3, -30, 30)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 5, 0.1, 100),
        cauchyprior(1, 1, -200, 200),
        lnormprior(log(1.6), 0.4214036, 0, 18),
        lnormprior(0.2, 0.5, 0.2, 18),
        normprior(0, 3, -30, 30)
      )
    }
    if(dmodel == 8){
      prior$priors[2,] <- lnormprior(1,2,0.2,20)
    }
  } else if (dmodel %in% c(19)){ #efsa family 1b, d not in exponent #EFSA-Gama
    if (dvariance == 1) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 2, 0.2, 20),
        cauchyprior(1, 1, -200, 200),
        lnormprior(0.2, 0.5, 0.2, 18),
        normprior(0, 3, -30, 30)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 2, 0.2, 20),
        cauchyprior(1, 1, -200, 200),
        lnormprior(0.2, 0.5, 0.2, 18),
        lnormprior(0, 2, 0, 100),
        normprior(-2, 3, -30, 30)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        cauchyprior(0, 1, -100, 100),
        lnormprior(1, 2, 0.2, 20),
        cauchyprior(1, 1, -200, 200),
        lnormprior(0.2, 0.5, 0.2, 18),
        normprior(0, 3, -30, 30)
      )
    }
  }

  return(prior)
}

##############################################################
# Standard Dichtomous
##############################################################
.bayesian_prior_dich <- function(model, degree = 2) {
  dmodel <- which(model == c(
    "hill", "gamma", "logistic", "log-logistic",
    "log-probit", "multistage", "probit",
    "qlinear", "weibull"
  ))

  if (dmodel == 1) { # HILL
    prior <- create_prior_list(
      normprior(-1, 2, -40, 40),
      normprior(0, 3, -40, 40),
      normprior(-3, 3.3, -40, 40),
      lnormprior(0.693147, 0.5, 0, 40)
    )
    prior <- create_dichotomous_prior(prior, "hill")
  }
  if (dmodel == 2) { # GAMMA
    prior <- create_prior_list(
      normprior(0, 2, -18, 18),
      lnormprior(0.693147180559945, 0.424264068711929, 0.2, 20),
      lnormprior(0, 1, 0, 1e4)
    )
    prior <- create_dichotomous_prior(prior, "gamma")
  }
  if (dmodel == 3) { # LOGISTIC
    prior <- create_prior_list(
      normprior(0, 1, -20, 20),
      lnormprior(0, 2, 0, 40)
    )
    prior <- create_dichotomous_prior(prior, "logistic")
  }
  if (dmodel == 4) { # LOG-LOGISTIC
    prior <- create_prior_list(
      normprior(0, 2, -20, 20),
      normprior(0, 1, -40, 40),
      lnormprior(0.693147180559945, 0.5, 0, 20)
    )
    prior <- create_dichotomous_prior(prior, "log-logistic")
  }
  if (dmodel == 5) { # LOG-PROBIT
    prior <- create_prior_list(
      normprior(0, 2, -20, 20),
      normprior(0, 1, -40, 40),
      lnormprior(0.693147180559945, 0.5, 0, 20)
    )
    prior <- create_dichotomous_prior(prior,"log-probit")
  }

  if (dmodel == 6) { # MULTISTAGE
    startP <- create_prior_list(
      normprior(0, 2, -20, 20),
      lnormprior(0, 0.5, 0, 100)
    )
    degree <- floor(degree)
    if (degree >= 2) { # make sure it is a positive degree
      for (ii in (2:degree)) {
        startP <- .combine_prior_lists(startP, lnormprior(0, 1, 0, 1e6))
      }
    }
    prior <- startP
    prior <- create_dichotomous_prior(prior, "multistage")
  }
  if (dmodel == 7) { # PROBIT
    prior <- create_prior_list(
      normprior(0, 1, -20, 20),
      lnormprior(0, 2, 0, 40)
    )
    prior <- create_dichotomous_prior(prior, "probit")
  }
  if (dmodel == 8) { # QLINEAR
    prior <- create_prior_list(
      normprior(0, 2, -20, 20),
      lnormprior(0, 1, 0, 18)
    )
    prior <- create_dichotomous_prior(prior, "qlinear")
  }
  if (dmodel == 9) { # WEIBULL
    prior <- create_prior_list(
      normprior(0, 2, -20, 20),
      lnormprior(0.424264068711929, 0.5, 0, 40),
      lnormprior(0, 1.5, 0, 1e4)
    )
    prior <- create_dichotomous_prior(prior, "weibull")
  }

  return(prior)
}

################################################# 33
# bayesian_prior_dich(model,variance)
##################################################
.MLE_bounds_continuous <- function(model, variance, degree = 2, is_increasing) {
  # dmodel <- which(model == c("hill", "exp-3", "exp-5", "power", "polynomial"))
  dmodel <- which(model == .continuous_models[-5]) #remove FUNL
  dvariance <- which(variance == c("normal", "normal-ncv", "lognormal"))

  # POLYNOMIAL *BLAH*
  if (dmodel == 5) {
    if (dvariance == 1) {
      prior <- create_prior_list(normprior(0, 5, -100000, 100000))

      for (ii in 1:degree) {
        if (is_increasing) {
          if (ii <= 2) {
            prior <- .combine_prior_lists(prior, normprior(0, 5, 0, 1e6))
          } else {
            prior <- .combine_prior_lists(prior, normprior(0, 5, -1e6, 1e6))
          }
        } else {
          if (ii <= 2) {
            prior <- .combine_prior_lists(prior, normprior(0, 5, -1e6, 0))
          } else {
            prior <- .combine_prior_lists(prior, normprior(0, 5, -1e6, 1e6))
          }
        }
      }

      prior <- .combine_prior_lists(prior, create_prior_list(normprior(0, 1, -18, 18)))
      prior[[1]][, 1] <- 0
      return(prior)
    }

    if (dvariance == 2) {
      prior <- create_prior_list(normprior(0, 5, 0, 1000))

      for (ii in 1:degree) {
        if (is_increasing) {
          if (ii <= 2) {
            prior <- .combine_prior_lists(prior, normprior(0, 5, 0, 1e6))
          } else {
            prior <- .combine_prior_lists(prior, normprior(0, 5, -1e6, 1e6))
          }
        } else {
          if (ii <= 2) {
            prior <- .combine_prior_lists(prior, normprior(0, 5, -1e6, 0))
          } else {
            prior <- .combine_prior_lists(prior, normprior(0, 5, -1e6, 1e6))
          }
        }
      }
      prior <- .combine_prior_lists(
        prior,
        create_prior_list(
          normprior(0, 1, 0, 18),
          normprior(0, 1, -18, 18)
        )
      )
      prior[[1]][, 1] <- 0
      return(prior)
    }
    if (dvariance == 3) {
      stop("Polynomial-Log-normal models are not allowed. Please
            choose normal or normal non-constant variance.\n")
    }
    return(prior)
  }

  # Hill
  if (dmodel == 1) {
    if (dvariance == 1) { # normal
      prior <- create_prior_list(normprior(0, 1, -100, 100))
      if (is_increasing) {
        prior <- .combine_prior_lists(prior, normprior(0, 2, 0, 100))
      } else {
        prior <- .combine_prior_lists(prior, normprior(0, 2, -100, 0))
      }
      prior <- .combine_prior_lists(prior, lnormprior(0, 1, 0, 5))
      prior <- .combine_prior_lists(prior, lnormprior(1, 1.2, 1, 18))
      prior <- .combine_prior_lists(prior, normprior(0, 1, -18, 18))
    } else if (dvariance == 2) { # normal ncv
      prior <- create_prior_list(normprior(0, 1, -100, 100))
      if (is_increasing) {
        prior <- .combine_prior_lists(prior, normprior(0, 2, 0, 100))
      } else {
        prior <- .combine_prior_lists(prior, normprior(0, 2, -100, 0))
      }
      prior <- .combine_prior_lists(prior, normprior(0, 1, 0, 5))
      prior <- .combine_prior_lists(prior, lnormprior(log(1.2), 1, 1, 18))
      prior <- .combine_prior_lists(prior, normprior(0, 2, -18, 18))
      prior <- .combine_prior_lists(prior, normprior(0, 2, -18, 18))
    } else if (dvariance == 3) { # log normal
      prior <- create_prior_list(normprior(0, 1, -100, 100))
      if (is_increasing) {
        prior <- .combine_prior_lists(prior, normprior(0, 2, 0, 100))
      } else {
        prior <- .combine_prior_lists(prior, normprior(0, 2, -100, 0))
      }
      prior <- .combine_prior_lists(prior, lnormprior(0, 1, 0, 5))
      prior <- .combine_prior_lists(prior, lnormprior(0, 1, 0, 18))
      prior <- .combine_prior_lists(prior, normprior(0, 1, -18, 18))
    }
  }

  # Exponential-3
  if (dmodel == 2) {
    if (dvariance == 2) {
      prior <- create_prior_list(
        normprior(0, 1, -0, 100),
        lnormprior(0, 0.5, 0, 100),
        normprior(0, 0.5, -20, 20), # log(c)
        lnormprior(0, 0.3, 1, 18), # d
        lnormprior(0, 0.5, 0, 18),
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 1) {
      prior <- create_prior_list(
        normprior(0, 0.1, 0, 100), # a
        lnormprior(0, 1, 0, 100), # b
        normprior(0, 0.5, -20, 20), # log(c)
        lnormprior(1, 0.2, 1, 18), # d
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        normprior(0, 0.1, -1000, 1000), # a
        lnormprior(0, 1, 0, 100), # b
        normprior(0, 1, -20, 20), # log(c)
        lnormprior(1, 0.2, 1, 18), # d
        normprior(0, 1, -18, 18)
      )
    }
  }

  # Exponential-5
  if (dmodel == 3) {
    if (dvariance == 1) {
      prior <- create_prior_list(
        normprior(0, 0.1, -100, 100), # a
        lnormprior(0, 1, 0, 100), # b
        normprior(0, 1, -10, 10), # log(c)
        lnormprior(1, 0.2, 0, 18), # d
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        normprior(0, 0.1, -100, 100),
        normprior(0, 1, 0, 100),
        normprior(0, 0.5, -20, 20), # log(c)
        lnormprior(0, 0.2, 1, 18), # d
        lnormprior(0, 0.5, 0, 18),
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 3) {
      prior <- create_prior_list(
        normprior(0, 0.1, -100, 100), # a
        lnormprior(0, 1, 0, 100), # b
        normprior(0, 1, -10, 10), # log(c)
        lnormprior(1, 0.2, 1, 18), # d
        normprior(0, 1, -18, 18)
      )
    }
  }

  # Power
  if (dmodel == 4) {
    if (dvariance == 1) {
      prior <- create_prior_list(
        normprior(0, 0.1, -100, 100), # a
        normprior(0, 1, -1e2, 1e2), # b
        lnormprior(1, 0.2, 0, 18), # k
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 2) {
      prior <- create_prior_list(
        normprior(0, 1, -100, 100), # a
        normprior(0, 1, -1e4, 1e4), # b
        lnormprior(0, 0.2, 1, 18), # k
        lnormprior(0, 0.250099980007996, 0, 18),
        normprior(0, 1, -18, 18)
      )
    } else if (dvariance == 3) {
      stop("Power-Log-normal models are not allowed. Please choose normal or normal non-constant variance. \n")
    }
  }
  
  if (dmodel %in% c(6:7, 10, 13,18)){ #4 parameter Aerts models
    prior <- create_prior_list(cauchyprior(0, 10, -1000, 1000),
                               lnormprior(0, 3, 1e-6, 100))
    if(is_increasing){
      prior <- .combine_prior_lists(prior, cauchyprior(1, 10, 0, 200))
    }else{
      prior <- .combine_prior_lists(prior, cauchyprior(-1, 10, -200, 0))
    }
    
    prior <- .combine_prior_lists(prior, lnormprior(log(1.6), 0.4214036, 0, 18))
    
    if (dvariance == 1) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    } else if (dvariance == 2) {
      prior <- .combine_prior_lists(prior, lnormprior(0, 2, 0, 100))
      prior <- .combine_prior_lists(prior, normprior(-2, 3, -30, 30))
    } else if (dvariance == 3) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    }
  } 
  if (dmodel %in% c(16,17)){ #4 parameter Aerts models (probit/logistic)
    prior <- create_prior_list(cauchyprior(0, 10, -1000, 1000),
                               normprior(0, 3, -100, 100))
    
    prior <- .combine_prior_lists(prior, cauchyprior(0, 10, -200, 200))
    prior <- .combine_prior_lists(prior, lnormprior(log(1.6), 0.4214036, 0, 18))
    
    if (dvariance == 1) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    } else if (dvariance == 2) {
      prior <- .combine_prior_lists(prior, lnormprior(0, 2, 0, 100))
      prior <- .combine_prior_lists(prior, normprior(-2, 3, -30, 30))
    } else if (dvariance == 3) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    }
  } 
  if (dmodel %in% c(8:9,11:12,14:15)){ #5 parameter Aerts models
    prior <- create_prior_list(cauchyprior(0, 10, -1000, 1000),
                               lnormprior(1, 5, 0.1, 100))
    if(dmodel == 8){ #gamma has numerical issues for large B
      prior <- create_prior_list(cauchyprior(0, 10, -1000, 1000),
                                 lnormprior(1, 5, 0.2, 20))
    }
    # if(dvariance == 3){#a has to be positive for lognormal models!
    #   prior <- create_prior_list(lnormprior(1, 10, 1e-6, 100),
    #                              lnormprior(1, 2, 0.2, 20))
    # }
    if(is_increasing){
      prior <- .combine_prior_lists(prior, cauchyprior(1, 5, 0, 200))
    }else{
      prior <- .combine_prior_lists(prior, cauchyprior(-1, 5, -200, 0))
    }
    
    prior <- .combine_prior_lists(prior, lnormprior(log(1.6), 0.4214036, 0, 18))
    prior <- .combine_prior_lists(prior, lnormprior(0.2, 0.5, 1e-6, 18))
    
    if (dvariance == 1) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    } else if (dvariance == 2) {
      prior <- .combine_prior_lists(prior, lnormprior(0, 2, 0, 100))
      prior <- .combine_prior_lists(prior, normprior(-2, 3, -30, 30))
    } else if (dvariance == 3) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    }
    
    # if(dmodel == 8){ #gamma has numerical issues for large B
    #   prior <- create_prior_list(normprior(1, 1, -100, 100),
    #                              lnormprior(0, 2, 0, 15),
    #                              normprior(1, 2, -20, 20),
    #                              lnormprior(log(1.6), 0.4214036, 0, 10),
    #                              lnormprior(0.2, 0.5, 0.1, 5),
    #                              normprior(0, 1, -10, 10)
    #                              )
    # }
  }
  if(dmodel %in% c(19)){ #efsa family 1b, d not in exponent

    #gamma has numerical issues for large B
    prior <- create_prior_list(cauchyprior(0, 10, -1000, 1000),
                               lnormprior(1, 5, 0.2, 20))
    if(is_increasing){
      prior <- .combine_prior_lists(prior, cauchyprior(1, 5, 0, 200))
    }else{
      prior <- .combine_prior_lists(prior, cauchyprior(-1, 5, -200, 0))
    }
    
    prior <- .combine_prior_lists(prior, lnormprior(0.2, 0.5, 1e-6, 18))
    
    if (dvariance == 1) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    } else if (dvariance == 2) {
      prior <- .combine_prior_lists(prior, lnormprior(0, 2, 0, 100))
      prior <- .combine_prior_lists(prior, normprior(-2, 3, -30, 30))
    } else if (dvariance == 3) {
      prior <- .combine_prior_lists(prior, normprior(0, 3, -30, 30))
    }
  }
  # if (dmodel %in% c(6:7, 10, 13, 16:17)){ #4 parameter Aerts models
  #   if (dvariance == 1) {
  #     prior <- create_prior_list(
  #       normprior(1, 1, -100, 100),
  #       lnormprior(0, 2, 0, 100),
  #       normprior(0, 2, -200, 200),
  #       lnormprior(log(1.6), 0.4214036, 0, 18),
  #       normprior(0, 1, -30, 30)
  #     )
  #   } else if (dvariance == 2) {
  #     prior <- create_prior_list(
  #       normprior(1, 1, -100, 100),
  #       lnormprior(0, 2, 0, 100),
  #       normprior(0, 2, -200, 200),
  #       lnormprior(log(1.6), 0.4214036, 0, 18),
  #       lnormprior(0, 1, 0, 100),
  #       normprior(-3, 1, -18, 18)
  #     )
  #   } else if (dvariance == 3) {
  #     prior <- create_prior_list(
  #       normprior(1, 1, -100, 100),
  #       lnormprior(0, 2, 0, 100),
  #       normprior(0, 2, -200, 200),
  #       lnormprior(log(1.6), 0.4214036, 0, 18),
  #       normprior(0, 1, -18, 18)
  #     )
  #   }
  # } else if (dmodel %in% c(8:9,11:12,14:15)){ #5 parameter Aerts models
  #   if (dvariance == 1) {
  #     prior <- create_prior_list(
  #       normprior(1, 1, -100, 100),
  #       lnormprior(0, 2, 0, 100),
  #       normprior(0, 2, -200, 200),
  #       lnormprior(log(1.6), 0.4214036, 0, 18),
  #       lnormprior(0.2, 0.5, 0, 18),
  #       normprior(0, 1, -30, 30)
  #     )
  #   } else if (dvariance == 2) {
  #     prior <- create_prior_list(
  #       normprior(1, 1, -100, 100),
  #       lnormprior(0, 2, 0, 100),
  #       normprior(0, 2, -200, 200),
  #       lnormprior(log(1.6), 0.4214036, 0, 18),
  #       lnormprior(0.2, 0.5, 0, 18),
  #       lnormprior(0, 1, 0, 100),
  #       normprior(-3, 1, -18, 18)
  #     )
  #   } else if (dvariance == 3) {
  #     prior <- create_prior_list(
  #       normprior(1, 1, -100, 100),
  #       lnormprior(0, 2, 0, 100),
  #       normprior(0, 2, -200, 200),
  #       lnormprior(log(1.6), 0.4214036, 0, 18),
  #       lnormprior(0.2, 0.5, 0, 18),
  #       normprior(0, 1, -18, 18)
  #     )
  #   }
  # }
  
  

  prior <- unclass(prior)
  prior[[1]][, 1] <- 0
  return(prior)
}
