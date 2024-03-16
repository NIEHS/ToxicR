
#' Fit a single dichotomous dose-response model to data.
#'
#' @param D A numeric vector or matrix of doses.
#' @param Y A numeric vector or matrix of responses.
#' @param N A numeric vector or matrix of the number of replicates at a dose.
#' @param model_type The mean model for the dichotomous model fit.  It can be one of the following: \cr
#'    "hill","gamma","logistic", "log-logistic", "log-probit"  ,"multistage"  ,"probit","qlinear","weibull"
#' @param fit_type the method used to fit (laplace, mle, or mcmc)
#' @param prior Used if you want to specify a prior for the data.
#' @param BMR This option specifies the benchmark response BMR. The BMR is defined in relation to the BMD calculation requested (see BMD).  By default, the "BMR = 0.1."
#' @param alpha Alpha is the specified nominal coverage rate for computation of the lower bound on the BMDL and BMDU, i.e., one computes a \eqn{100\times(1-\alpha)\%} .  For the interval (BMDL,BMDU) this is a \eqn{100\times(1-2\alpha)\% } confidence interval.  By default, it is set to 0.05.
#' @param degree the number of degrees of a polynomial model. Only used for polynomial models.
#' @param samples the number of samples to take (MCMC only)
#' @param burnin the number of burnin samples to take (MCMC only)
#' @param threads specify the number of OpenMP threads to use for the calculations. Default = 2
#' @param seed specify the GSL seed. Default = 12331
#'
#' @return Returns a model object class with the following structure:
#' \itemize{
#'    \item \code{full_model}:  The model along with the likelihood distribution.
#'    \item \code{parameters}: The parameter estimates produced by the procedure, which are relative to the model '
#'                             given in \code{full_model}.  The last parameter is always the estimate for \eqn{\log(\sigma^2)}.
#'    \item \code{covariance}: The variance-covariance matrix for the parameters.
#'    \item \code{bmd_dist}:  Quantiles for the BMD distribution.
#'    \item \code{bmd}:  A vector containing the benchmark dose (BMD) and \eqn{100\times(1-2\alpha)} confidence intervals.
#'    \item \code{maximum}:  The maximum value of the likelihod/posterior.
#'    \item \code{gof_p_value}:  GOF p-value for the Pearson \eqn{\chi^2} GOF test.
#'    \item \code{gof_chi_sqr_statistic}: The GOF statistic.
#'    \item \code{prior}:     This value gives the prior for the Bayesian analysis.
#'    \item \code{model}:     Parameter specifies t mean model used.
#'    \item \code{data}:      The data used in the fit.
#'    \itemize{
#'        When MCMC is specified, an additional variable \code{mcmc_result}
#'        has the following two variables:
#'        \item \code{PARM_samples}:  matrix of parameter samples.
#'        \item \code{BMD_samples}: vector of BMD sampled values.
#'    }
#' }
#'
#' @examples
#' mData <- matrix(c(
#'   0, 2, 50,
#'   1, 2, 50,
#'   3, 10, 50,
#'   16, 18, 50,
#'   32, 18, 50,
#'   33, 17, 50
#' ), nrow = 6, ncol = 3, byrow = TRUE)
#' D <- mData[, 1]
#' Y <- mData[, 2]
#' N <- mData[, 3]
#' model <- single_dichotomous_fit(D, Y, N, model_type = "hill", fit_type = "laplace")
#' summary(model)
#'
single_dichotomous_fit <- function(D, Y, N, model_type, fit_type = "laplace",
                                   prior = NULL, BMR = 0.1,
                                   alpha = 0.05, degree = 2, samples = 21000,
                                   burnin = 1000, threads=2, seed = 12331) {
  .setseedGSL(seed)
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  N <- as.matrix(N)

  DATA <- cbind(D, Y, N)
  test <- .check_for_na(DATA)
  Y <- Y[test == TRUE, , drop = F]
  D <- D[test == TRUE, , drop = F]
  N <- N[test == TRUE, , drop = F]

  if (is.null(prior)) {
    prior <- .bayesian_prior_dich(model_type, degree)
  } else {
    if (!("BMD_Bayes_dichotomous_model" %in% class(prior))) {
      stop("Prior is not correctly specified.")
    }
    model_type <- prior$mean
    if (model_type == "multistage") {
      degree <- prior$degree
    }
    if (fit_type == "mle") {
      stop("A Bayesian prior model was specified, but MLE was requested.")
    }
  }

  dmodel <- which(model_type == c(
    "hill", "gamma", "logistic", "log-logistic",
    "log-probit", "multistage", "probit",
    "qlinear", "weibull"
  ))
  DATA <- cbind(D, Y, N)
  o1 <- c(BMR, alpha, -9999)
  o2 <- c(1, degree)

  if (identical(dmodel, integer(0))) {
    stop('Please specify one of the following model types:
            "hill", "gamma", "logistic", "log-logistic"
            "log-probit", "multistage"
            "probit", "qlinear", "weibull"')
  }
  if (dmodel == 6) {
    if ((o2[2] < 2) + (o2[2] > nrow(DATA) - 1) > 0) {
      stop("The multistage model needs to have between
               2 and nrow(DATA)-1 paramaters. If degree = 1
               use the quantal linear model.")
    }
  }
  fit_type = tolower(fit_type)
  fitter <- which(fit_type == c("mle", "laplace", "mcmc"))
  if (identical(fitter, integer(0))) {
    stop('The fit_type variable must be either "laplace","mle", or "mcmc"\n')
  }


  if (fitter == 1) { # MLE fit
    bounds <- .bmd_default_frequentist_settings(model_type, degree)
    .set_threads(threads)
    temp <- .run_single_dichotomous(dmodel, DATA, bounds, o1, o2)
    # class(temp$bmd_dist) <- "BMD_CDF"
    temp_me <- temp$bmd_dist

    temp_me <- temp_me[!is.infinite(temp_me[, 1]), ]
    temp_me <- temp_me[!is.na(temp_me[, 1]), ]
    temp_me <- temp_me[!is.nan(temp_me[, 1]), ]
    if (is.null(nrow(temp_me))){temp$bmd <- c(temp$bmd, NA, NA)} else if(nrow(temp_me) > 5) {
      te <- splinefun(temp_me[, 2], temp_me[, 1], method = "monoH.FC",ties=mean)
      temp$bmd <- c(temp$bmd, te(alpha), te(1 - alpha))
    } else {
      temp$bmd <- c(temp$bmd, NA, NA)
    }
    temp$bounds <- bounds
    temp$model <- model_type
    temp$data <- DATA
    temp$options <- c(BMR, alpha, samples, burnin)
    class(temp) <- "BMDdich_fit_maximized"
  }

  if (fitter == 2) { # laplace fit
    .set_threads(threads)
    temp <- .run_single_dichotomous(dmodel, DATA, prior$priors, o1, o2)
    # class(temp$bmd_dist) <- "BMD_CDF"
    temp_me <- temp$bmd_dist
    temp_me <- temp_me[!is.infinite(temp_me[, 1]), ]
    temp_me <- temp_me[!is.na(temp_me[, 1]), ]
    temp_me <- temp_me[!is.nan(temp_me[, 1]), ]
    if (is.null(nrow(temp_me))){temp$bmd <- c(temp$bmd, NA, NA)} else if(nrow(temp_me) > 5) {
      te <- splinefun(temp_me[, 2], temp_me[, 1], method = "monoH.FC",ties=mean)
      temp$bmd <- c(temp$bmd, te(alpha), te(1 - alpha))
    } else {
      temp$bmd <- c(temp$bmd, NA, NA)
    }
    temp$prior <- prior
    temp$model <- model_type
    temp$data <- DATA
    temp$options <- c(BMR, alpha, samples, burnin)
    class(temp) <- "BMDdich_fit_maximized"
  }
  if (fitter == 3) {
    .set_threads(threads)
    temp <- .run_dichotomous_single_mcmc(
      dmodel, DATA[, 2:3, drop = F], DATA[, 1, drop = F], prior$priors,
      c(BMR, alpha, samples, burnin)
    )
    # class(temp$fitted_model$bmd_dist) <- "BMD_CDF"
    temp$bmd_dist <- cbind(quantile(temp$mcmc_result$BMD_samples, seq(0.005, 0.995, 0.005),na.rm=TRUE), seq(0.005, 0.995, 0.005))

    temp$options <- options <- c(BMR, alpha, samples, burnin)
    temp$prior <- prior <- list(prior = prior)
    temp$model <- model_type
    temp$data <- DATA
    temp$full_model <- temp$fitted_model$full_model
    temp$parameters <- temp$fitted_model$parameters
    temp$covariance <- temp$fitted_model$covariance
    temp$maximum <- temp$fitted_model$maximum
    temp$bmd <- as.numeric(c(mean(temp$mcmc_result$BMD_samples,na.rm=TRUE), quantile(temp$mcmc_result$BMD_samples, c(alpha, 1 - alpha), na.rm = TRUE)))
    temp$fitted_model <- NULL
    class(temp) <- "BMDdich_fit_MCMC"
  }
  names(temp$bmd) <- c("BMD", "BMDL", "BMDU")
  return(temp)
}




.bmd_default_frequentist_settings <- function(model, degree = 2) {
  dmodel <- which(model == c(
    "hill", "gamma", "logistic", "log-logistic",
    "log-probit", "multistage", "probit",
    "qlinear", "weibull"
  ))
  if (dmodel == 1) { # HILL
    prior <- matrix(c(
      0, 0, 2, -18, 18,
      0, 0, 2, -18, 18,
      0, 0, 0.5, -18, 18,
      0, 1, 0.250099980007996, 1.00E+00, 18
    ), nrow = 4, ncol = 5, byrow = T)
  }
  if (dmodel == 2) { # GAMMA
    prior <- matrix(c(
      0, 0, 2, -18, 18,
      0, 1.5, 0.424264068711929, 1, 18,
      0, 1, 1, 0, 1000
    ), nrow = 3, ncol = 5, byrow = T)
  }
  if (dmodel == 3) { # LOGISTIC
    prior <- matrix(c(
      0, -2, 2, -18, 18,
      0, 0.1, 1, 0.00E+00, 1e4
    ), nrow = 2, ncol = 5, byrow = T)
  }
  if (dmodel == 4) { # LOG-LOGISTIC
    prior <- matrix(c(
      0, -1.65, 2, -18, 18,
      0, -2, 1, -18, 18,
      0, 1.2, 0.5, 0, 1e4
    ), nrow = 3, ncol = 5, byrow = T)
  }
  if (dmodel == 5) { # LOG-PROBIT
    prior <- matrix(c(
      0, 0.01, 2, -18, 18,
      0, 0, 1, -8, 8,
      0, 1, 0.5, 0, 1000
    ), nrow = 3, ncol = 5, byrow = T)
  }

  if (dmodel == 6) { # MULTISTAGE
    temp <- matrix(c(
      0, -2, 2, -18, 18,
      0, 1, 0.25, 0, 1e4,
      0, 0.1, 1, 0, 1.00E+06
    ), nrow = 3, ncol = 5, byrow = T)
    prior <- matrix(c(0, 0.1, 1, 0, 1.00E+06), nrow = 1 + degree, ncol = 5, byrow = T)
    prior[1:3, ] <- temp
  }
  if (dmodel == 7) { # PROBIT
    prior <- matrix(c(
      0, -2, 2, -8, 8,
      0, 0.1, 1, 0.00E+00, 1000
    ), nrow = 2, ncol = 5, byrow = T)
  }
  if (dmodel == 8) { # QLINEAR
    prior <- matrix(c(
      0, -2, 2, -18, 18,
      0, 1, 1, 0, 1000
    ), nrow = 2, ncol = 5, byrow = T)
  }
  if (dmodel == 9) { # WEIBULL
    prior <- matrix(c(
      0, 0, 2, -18, 18,
      0, 1, 1, 0, 50,
      0, 1, 0.424264068711929, 1.00E-06, 1000
    ), nrow = 3, ncol = 5, byrow = T)
  }

  return(prior)
}
