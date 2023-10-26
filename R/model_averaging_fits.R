#' Fit a model averaged continuous BMD model.
#'
#' @title ma_continuous_fit - Fit a model averaged continuous BMD model.
#' @param D doses matrix
#' @param Y response matrix
#' @param model_list a list of configurations for the single models (priors, model type).  To create a model list, one creates a list of
#'                   continuous model priors using \code{create_continuous_prior}.
#' @param fit_type the method used to fit ("laplace", "mle", or "mcmc")
#' @param BMR_TYPE Specifies the type of benchmark dose analysis to be performed. For continuous models, there are four types of BMD definitions that are commonly used:
#'    \itemize{
#'      \item Standard deviation is the default option, but it can be explicitly specified with 'BMR_TYPE = "sd"' This definition defines the BMD as the dose associated with the mean/median changing a specified number of standard deviations from the mean at the control dose., i.e., it is the dose, BMD, that solves \eqn{\mid f(dose)-f(0) \mid = BMR \times \sigma}
#'      \item Relative deviation can be specified with 'BMR_TYPE = "rel"'. This defines the BMD as the dose that changes the control mean/median a certain percentage from the background dose, i.e. it is the dose, BMD that solves \eqn{\mid f(dose) - f(0) \mid = (1 \pm BMR) f(0)} 
#'      \item Hybrid deviation can be specified with 'BMR_TYPE = "hybrid"'.  This defines the BMD that changes the probability of an adverse event by a stated amount relitive to no exposure (i.e 0).  That is, it is the dose, BMD, that solves \eqn{\frac{Pr(X > x| dose) - Pr(X >x|0)}{Pr(X < x|0)} = BMR}. For this definition, \eqn{Pr(X < x|0) = 1 - Pr(X > X|0) = \pi_0}, where \eqn{0 \leq \pi_0 < 1} is defined by the user as "point_p," and it defaults to 0.01.  Note: this discussion assumed increasing data.  The fitter determines the direction of the data and inverts the probability statements for decreasing data.
#'      \item Absolute deviation can be specified with 'BMR_TYPE="abs"'. This defines the BMD as an absolute change from the control dose of zero by a specified amount. That is the BMD is the dose that solves the equation \eqn{\mid f(dose) - f(0) \mid = BMR}. 
#'    }  
#' @param BMR This option specifies the benchmark response BMR. The BMR is defined in relation to the BMD calculation requested (see BMD).  By default, the "BMR = 0.1."
#' @param point_p This option is only used for hybrid BMD calculations. It defines a probability that is the cutpoint for observations.  It is the probability that observations have this probability, or less, of being observed at the background dose.
#' @param alpha Alpha is the specified nominal coverage rate for computation of the lower bound on the BMDL and BMDU, i.e., one computes a \eqn{100\times(1-\alpha)\%} confidence interval.  For the interval (BMDL,BMDU) this is a \eqn{100\times(1-2\alpha)\% }.  By default, it is set to 0.05.
#' @param EFSA This option is a logical allowing users to use the default EFSA models to average instead of the original option in previous ToxicR versions. If false, the models are: 
#'        hill, hill normal-ncv, exp-3, exp-3 normal-ncv, exp-3 lognormal, exp-5, exp-5 normal-ncv, exp-5 lognormal, power, power normal-ncv. If true, it follows the EFSA guidelines;
#'        these models are (normal then lognormal for each model): exp-aerts, invexp-aerts, hill-aerts, lognormal-aerts, gamma-efsa, LMS, probit-aerts, and logistic-aerts. 
#'        See \code{\link{single_continuous_fit}} for details on the models.
#' @param samples the number of samples to take (MCMC only)
#' @param burnin the number of burnin samples to take (MCMC only)
#' @param BMD_TYPE Deprecated version of BMR_TYPE that specifies the type of benchmark dose analysis to be performed
#' @param threads specify the number of OpenMP threads to use for the calculations. Default = 2
#' @return This function model object containing a list of individual fits and model averaging fits
#' \itemize{
#'  \item \code{Individual_Model_X}: Here \code{X} is a number \eqn{1\leq X \leq n,} where \eqn{n}
#'         is the number of models in the model average.  For each \code{X}, this is an individual model
#'         fit identical to what is returned in `\code{single_continuous_fit}.'
#'  \item \code{ma_bmd}: The CDF of the model averaged BMD distribution.
#'  \item \code{posterior_probs}: The posterior model probabilities used in the MA.
#'  \item \code{bmd}: The BMD and the \eqn{100\times(1-2\alpha)\%} confidence intervals.
#' }
#' @examples
#' \donttest{
#' hill_m <- function(doses) {
#'   returnV <- 481 - 250.3 * doses^1.3 / (40^1.3 + doses^1.3)
#'   return(returnV)
#' }
#' doses <- rep(c(0, 6.25, 12.5, 25, 50, 100), each = 10)
#' mean <- hill_m(doses)
#' y <- rnorm(length(mean), mean, 20.14)
#' model <- ma_continuous_fit(doses, y, fit_type = "laplace", BMR_TYPE = "sd", BMR = 1)
#' summary(model)
#' }
ma_continuous_fit <- function(D, Y, model_list = NA, fit_type = "laplace",
                              BMR_TYPE = "sd", BMR = 0.1, point_p = 0.01,
                              alpha = 0.05, EFSA = TRUE, samples = 21000,
                              burnin = 1000, BMD_TYPE = NA, threads=2) {
  myD <- Y
  Y <- as.matrix(Y)
  D <- as.matrix(D)

  is_neg <- .check_negative_response(Y)

  DATA <- cbind(D, Y)
  test <- .check_for_na(DATA)
  Y <- Y[test == TRUE, , drop = F]
  D <- D[test == TRUE, , drop = F]
  DATA <- cbind(D, Y)

  current_models <- c("hill", "exp-3", "exp-5", "power", "FUNL", "exp-aerts", "invexp-aerts", 
        "gamma-aerts", "invgamma-aerts", "hill-aerts", "lomax-aerts", "invlomax-aerts", "lognormal-aerts",
        "logskew-aerts", "invlogskew-aerts", "logistic-aerts", "probit-aerts", "LMS", "gamma-efsa")
  current_dists <- c("normal", "normal-ncv", "lognormal")
  fit_type = tolower(fit_type)
  type_of_fit <- which(fit_type == c("laplace", "mcmc"))

  if(!is.na(BMD_TYPE)){
    warning("BMD_TYPE is deprecated. Please use BMR_TYPE instead")
  }else{
    BMD_TYPE = BMR_TYPE
  }
  
  rt <- which(BMD_TYPE == c("abs", "sd", "rel", "hybrid"))
  if (rt == 4) {
    rt <- 6
  }
  if (identical(rt, integer(0))) {
    stop("Please specify one of the following BMRF types:
    		  'abs','sd','rel','hybrid'")
  }

  if (rt == 4) {
    rt <- 6
  } # internally hybrid is coded as 6


  if (is.na(model_list[[1]][1])) {
    # no prior distribution specified as a parameter
    if(!(EFSA)){
      model_list <- c(rep("hill", 2), rep("exp-3", 3), rep("exp-5", 3), rep("power", 2))
      distribution_list <- c(
        "normal", "normal-ncv", rep(c("normal", "normal-ncv", "lognormal"), 2),
        "normal", "normal-ncv"
      )
    }else{
      model_list <- c(rep("exp-aerts", 2), rep("invexp-aerts",2), rep("hill-aerts", 2), rep("lognormal-aerts",2),
      rep("gamma-efsa", 2), rep("LMS", 2), rep("probit-aerts", 2), rep("logistic-aerts",2))
      distribution_list <- rep(c("normal", "lognormal"), 8)
    }

    if (is_neg) {
      tmpIdx <- which(distribution_list == "lognormal")
      model_list <- model_list[-tmpIdx]
      distribution_list <- distribution_list[-tmpIdx]
      if (length(distribution_list) > 1) # need at least 2 models for model averaging
        {
          warning("Negative response values were found in the data.  All lognormal
        models were removed from the analysis.")
        } else {
        stop("Negative response values were found in the data.  All lognormal models were removed from the analysis, but there were not enough models available for the MA.")
      }
    }
    prior_list <- list()
    for (ii in 1:length(model_list)) {
      PR <- .bayesian_prior_continuous_default(model_list[ii], distribution_list[ii], 2)
      # specify variance of last parameter to variance of response
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
      t_prior_result <- create_continuous_prior(PR, model_list[ii], distribution_list[ii], 2)
      PR <- t_prior_result$prior
      prior_list[[ii]] <- list(prior = PR, model_tye = model_list[ii], dist = distribution_list[ii])
    }
    model_list2 <- model_list
  } else {
    prior_list <- list()
    distribution_list <- c()
    model_list2 <- c()
    for (ii in 1:length(model_list)) {
      temp_prior <- model_list[[ii]]


      if (!("BMD_Bayes_continuous_model" %in% class(temp_prior))) {
        stop("Prior is not the correct form. Please use a Bayesian Continuous Prior Model.")
      }
      result <- .parse_prior(temp_prior)
      distribution <- result$distribution
      model_type <- result$model

      if (model_type == "polynomial") {
        stop("Polynomial models are not allowed in model averaging.")
      }
      a <- list(
        model = model_type, dist = distribution,
        prior = result$prior
      )
      prior_list[[ii]] <- a
      distribution_list <- c(distribution_list, distribution)
      model_list2 <- c(model_list2, model_type)
    }
  }

  models <- rep(0, length(prior_list))
  dlists <- rep(0, length(prior_list))
  priors <- list()
  # permuteMat <- cbind(c(1, 2, 3, 4, 5), c(6, 3, 5, 8, 10)) # c++ internal adjustment
  permuteMat = cbind(c(1,2,3,4,5,7:20),c(6,3,5,8,10, 11:24))
  for (ii in 1:length(prior_list)) {
    models[ii] <- permuteMat[which(prior_list[[ii]]$model == current_models), 2] # readjust for c++ internal
    priors[[ii]] <- prior_list[[ii]]$prior
    dlists[ii] <- which(prior_list[[ii]]$dist == current_dists)
  }


  ###################
  DATA <- cbind(D, Y)
  if (ncol(DATA) == 4) {
    colnames(DATA) <- c("Dose", "Resp", "N", "StDev")
    sstat <- T
  } else if (ncol(DATA) == 2) {
    colnames(DATA) <- c("Dose", "Resp")
    sstat <- F
  } else {
    stop("The data do not appear to be in the correct format.")
  }

  model_data <- list()
  model_data$X <- D
  model_data$SSTAT <- DATA

  if (sstat == T) {
    temp.fit <- lm(model_data$SSTAT[, 2] ~ model_data$X,
      weights = (1 / model_data$SSTAT[, 4]^2) * model_data$SSTAT[, 3]
    )
  } else {
    temp.fit <- lm(model_data$SSTAT[, 2] ~ model_data$X)
  }

  # Determine if there is an increasing or decreasing trend for BMD
  is_increasing <- F
  if (coefficients(temp.fit)[2] > 0) {
    is_increasing <- T
  }
  if (!is_increasing) {
    if (BMD_TYPE == "rel") {
      BMR <- 1 - BMR
    }
  }
  options <- c(rt, BMR, point_p, alpha, is_increasing, samples, burnin)
  
  .set_threads(threads)
  if (fit_type == "mcmc") {
    temp_r <- .run_continuous_ma_mcmc(
      priors, models, dlists, Y, D,
      options
    )
    tempn <- temp_r$ma_results

    tempm <- temp_r$mcmc_runs
    # clean up the run

    temp <- list()
    idx <- grep("Fitted_Model", names(tempn))


    jj <- 1
    for (ii in idx) {
      temp[[jj]] <- list()
      temp[[jj]]$mcmc_result <- tempm[[ii]]
      # remove the unecessary 'c' column from the exponential fit
      if ("exp-3" %in% model_list[jj]) {
        temp[[jj]]$mcmc_result$PARM_samples <- temp[[jj]]$mcmc_result$PARM_samples[, -3]
      }
      temp[[jj]]$full_model <- tempn[[ii]]$full_model
      temp[[jj]]$parameters <- tempn[[ii]]$parameters
      temp[[jj]]$covariance <- tempn[[ii]]$covariance
      temp[[jj]]$maximum <- tempn[[ii]]$maximum
      temp[[jj]]$bmd_dist <- tempn[[ii]]$bmd_dist

      temp[[jj]]$prior <- priors[[which(ii == idx)]]
      temp[[jj]]$data <- cbind(D, Y)
      temp[[jj]]$model <- prior_list[[jj]]$model # tolower(trimws(gsub("Model: ","",temp[[ii]]$full_model)))

      data_temp <- temp[[jj]]$bmd_dist
      data_temp <- data_temp[!is.infinite(data_temp[, 1]) & !is.na(data_temp[, 1]), ]
      temp[[jj]]$bmd <- c(NA, NA, NA)


      if (length(data_temp) > 0) {
        ii <- nrow(data_temp)

        temp[[jj]]$bmd_dist <- data_temp
        if (length(data_temp) > 10) {
          te <- splinefun(data_temp[, 2, drop = F], data_temp[, 1, drop = F], method = "monoH.FC",ties=mean)
          temp[[jj]]$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
        }
      }
      #add NAs if bad hessian or NaN BMD
      if(det(temp[[jj]]$covariance) < 0 || is.na(temp[[jj]]$bmd[1])){
        tempn$posterior_probs[jj] <- NA
      }
      names(temp[[jj]]$bmd) <- c("BMD", "BMDL", "BMDU")
      temp[[jj]]$options <- options
      class(temp[[jj]]) <- "BMDcont_fit_MCMC"
      jj <- jj + 1
    }
    # for (ii in idx_mcmc)
    # names(temp) <- sprintf("Individual_Model_%s", 1:length(idx)) 
    names(temp) <- sprintf("Indiv_%s_%s", model_list2, distribution_list)
    # print(tempn)
    temp$ma_bmd <- tempn$ma_bmd

    data_temp <- temp$ma_bmd
    data_temp <- data_temp[!is.infinite(data_temp[, 1]) & !is.na(data_temp[, 1]), ]
    #    data_temp = data_temp[!is.na(data_temp[,1]),]
    temp$bmd <- c(NA, NA, NA)

    if(min(data_temp[,2]) > alpha){
        data_temp <- rbind(c(0,0), data_temp)
        warning("BMDL may be inaccurate")
    }
    if (length(data_temp) > 0) {
      ii <- nrow(data_temp)

      temp$ma_bmd <- data_temp
      #tempn$posterior_probs[is.nan(tempn$posterior_probs)] <- 0
      if (length(data_temp) > 10 && (abs(sum(tempn$posterior_probs,na.rm=TRUE) - 1) <= 1e-8)) {
        te <- splinefun(data_temp[, 2, drop = F], data_temp[, 1, drop = F], method = "monoH.FC",ties=mean)
        temp$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
      } else {
        # error with the posterior probabilities
        temp$bmd <- c(Inf, 0, Inf)
      }
    }

    names(temp$bmd) <- c("BMD", "BMDL", "BMDU")
    temp$posterior_probs <- tempn$posterior_probs
    names(temp$posterior_probs) <- paste(model_list2, distribution_list, sep="_")
    class(temp) <- c("BMDcontinuous_MA", "BMDcontinuous_MA_mcmc")
    return(temp)
  } else {
    temp <- .run_continuous_ma_laplace(
      priors, models, dlists, Y, D,
      options
    )
    t_names <- names(temp)

    idx <- grep("Fitted_Model", t_names)
    jj <- 1
    for (ii in idx) {
      temp[[ii]]$prior <- priors[[which(ii == idx)]]
      temp[[ii]]$data <- cbind(D, Y)
      temp[[ii]]$model <- prior_list[[jj]]$model

      data_temp <- temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[, 1]), ]
      if (length(data_temp) > 0) {
        data_temp <- data_temp[!is.na(data_temp[, 1]), ]
        if (nrow(data_temp) > 6) {
          te <- splinefun(sort(data_temp[, 2, drop = F]), sort(data_temp[, 1, drop = F]), method = "monoH.FC",ties=mean)
          temp[[ii]]$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
          if (max(data_temp[, 2]) < 1 - alpha) {
            temp[[ii]]$bmd[3] <- 1e300
          }
        } else {
          temp[[ii]]$bmd <- c(NA, NA, NA)
        }
      }

      #add NAs if bad hessian or NaN BMD
      if(det(temp[[ii]]$covariance) < 0 || is.na(temp[[ii]]$bmd[1])){
        temp$posterior_probs[ii] <- NA
      }

      names(temp[[ii]]$bmd) <- c("BMD", "BMDL", "BMDU")
      temp[[ii]]$options <- options
      # names(temp)[ii] <- sprintf("Individual_Model_%s", ii)
      names(temp)[ii] <- sprintf("Indiv_%s_%s",model_list2[ii], distribution_list[ii])
      class(temp[[ii]]) <- "BMDcont_fit_maximized"
      jj <- jj + 1
    }

    temp_me <- temp$ma_bmd
    temp$bmd <- c(NA, NA, NA)
    if (!is.null(dim(temp_me))) {
      temp_me <- temp_me[!is.infinite(temp_me[, 1]), ]
      temp_me <- temp_me[!is.na(temp_me[, 1]), ]
      temp_me <- temp_me[!is.nan(temp_me[, 1]), ]
      #temp$posterior_probs[is.nan(temp$posterior_probs)] <- 0

      if(min(temp_me[,2]) > alpha){
        temp_me <- rbind(c(0,0), temp_me)
        warning("BMDL may be inaccurate")
      }
      if ((nrow(temp_me) > 10) && abs(sum(temp$posterior_probs, na.rm=TRUE) - 1) <= 1e-8) {
        te <- splinefun(sort(temp_me[, 2, drop = F]), sort(temp_me[, 1, drop = F]), method = "monoH.FC",ties=mean)
        temp$bmd <- c(te(0.5), te(alpha), te(1 - alpha))

        if (max(temp_me[, 2]) < 1 - alpha) {
          temp$bmd[3] <- 1e300
        }
      } else {
        temp$bmd <- c(Inf, 0, Inf)
      }
    }
    temp$ma_bmd <- temp_me
    names(temp$bmd) <- c("BMD", "BMDL", "BMDU")
    #temp$posterior_probs <- temp$posterior_probs
    names(temp$posterior_probs) <- paste(model_list2, distribution_list, sep="_")
    class(temp) <- c("BMDcontinuous_MA", "BMDcontinuous_MA_laplace")
    return(temp)
  }
}

.dichotomous_model_type <- function(model_name) {
  # based upon the following enum in the c++ file bmdStruct.h
  # enum dich_model {d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
  # d_logprobit=5, d_multistage=6,d_probit=7,
  # d_qlinear=8,d_weibull=9};
  result <- which(model_name == c(
    "hill", "gamma", "logistic", "log-logistic", "log-probit", "multistage", "probit",
    "qlinear", "weibull"
  ))

  if ((identical(result, integer(0)))) {
    stop("The model requested is not defined.")
  }

  return(result)
}

#' Fit a model averaged dichotomous BMD model.
#'
#' @title ma_dichotomous_fit - Fit a model averaged dichotomous BMD model.
#' @param D doses matrix
#' @param Y response matrix
#' @param N number of replicates matrix
#' @param model_list a list of configurations for the single models (priors, model type)
#' @param fit_type the method used to fit (laplace, mle, or mcmc)
#' @param BMR_TYPE Specifies the type of benchmark dose analysis to be performed. For dichotomous models, "extra" is the default BMR type. Other strings lead to "added" BMR calculations.
#' @param BMR This option specifies the benchmark response BMR. The BMR is defined in relation to the BMD calculation requested (see BMD).  By default, the "BMR = 0.1."
#' @param point_p This option is only used for hybrid BMD calculations. It defines a probability that is the cutpoint for observations.  It is the probability that observations have this probability, or less, of being observed at the background dose.
#' @param alpha Alpha is the specified nominal coverage rate for computation of the lower bound on the BMDL and BMDU, i.e., one computes a \eqn{100\times(1-\alpha)\% }.  For the interval (BMDL,BMDU) this is a \eqn{100\times(1-2\alpha)\% }.  By default, it is set to 0.05.
#' @param samples the number of samples to take (MCMC only)
#' @param burnin the number of burnin samples to take (MCMC only)
#' @param BMD_TYPE Deprecated version of BMR_TYPE that specifies the type of benchmark dose analysis to be performed
#' @return a model object containing a list of single models
#' \itemize{
#'  \item \code{Individual_Model_X}: Here \code{X} is a number \eqn{1\leq X \leq n,} where \eqn{n}
#'         is the number of models in the model average.  For each \code{X}, this is an individual model
#'         fit identical to what is returned in `\code{single_continuous_fit}.'
#'  \item \code{ma_bmd}: The CDF of the model averaged BMD distribution.
#'  \item \code{posterior_probs}: The posterior model probabilities used in the MA.
#'  \item \code{bmd}: The BMD and the \eqn{100\times(1-2\alpha)\%} confidence intervals.
#' }
#' @examples
#' \donttest{
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
#' model <- ma_dichotomous_fit(D, Y, N)
#' 
#' summary(model)
#' }
#' @export
ma_dichotomous_fit <- function(D, Y, N, model_list = integer(0), fit_type = "laplace",
                               BMR_TYPE = "extra",
                               BMR = 0.1, point_p = 0.01, alpha = 0.05, samples = 21000,
                               burnin = 1000, BMD_TYPE = NA, threads = 2) {
  D <- as.matrix(D)
  Y <- as.matrix(Y)
  N <- as.matrix(N)
  fit_type = tolower(fit_type)
  DATA <- cbind(D, Y, N)
  test <- .check_for_na(DATA)
  Y <- Y[test == TRUE, , drop = F]
  D <- D[test == TRUE, , drop = F]
  N <- N[test == TRUE, , drop = F]

  priors <- list()
  temp_prior_l <- list()
  tmodel_list <- list()
  if (length(model_list) < 1) {
    model_list <- .dichotomous_models
    model_i <- rep(0, length(model_list))
    for (ii in 1:length(model_list)) {
      temp_prior_l[[ii]] <- .bayesian_prior_dich(model_list[ii])
      priors[[ii]] <- temp_prior_l[[ii]]$priors
      model_i[ii] <- .dichotomous_model_type(model_list[ii])
    }
  } else {
    if (!("list" %in% class(model_list))) {
      stop("Please pass a list of priors.")
    }
    tmodel_list <- model_list
    model_list <- rep("", length(model_list))
    model_i <- rep(0, length(model_list))
    for (ii in 1:length(model_list)) {
      if (!("BMD_Bayes_dichotomous_model" %in% class(tmodel_list[[ii]]))) {
        stop("One of the specified models is not a 'BMD_Bayes_dichotomous_model.'")
      }
      temp_prior_l <- tmodel_list[[ii]]
      priors[[ii]] <- temp_prior_l$priors
      model_list[ii] <- temp_prior_l$mean
      model_i[ii] <- .dichotomous_model_type(model_list[ii])
    }
  }

  distribution_list <- rep("", length(model_list))
  # return(list(priors,model_list,model_i))

  model_p <- rep(1, length(model_list)) / length(model_list) # background prior is even
  o1 <- c(BMR, alpha)

  if(!is.na(BMD_TYPE)){
    warning("BMD_TYPE is deprecated. Please use BMR_TYPE instead")
  }else{
    BMD_TYPE = BMR_TYPE
  }

  if (BMD_TYPE == "extra") {
    BTYPE <- 1
  } else {
    BTYPE <- 2 # Added risk
  }

  o2 <- c(BTYPE, 2, samples, burnin)

  data <- as.matrix(cbind(D, Y, N))
  .set_threads(threads)
  if (fit_type == "laplace") {
    # Laplace Run
    temp <- .run_ma_dichotomous(
      data, priors, model_i,
      model_p, FALSE, o1, o2
    )
    # clean up the run
    temp$ma_bmd <- temp$BMD_CDF
    # TO DO : DELETE temp$BMD_CDF
    te <- splinefun(temp$ma_bmd[!is.infinite(temp$ma_bmd[, 1]), 2],
      temp$ma_bmd[!is.infinite(temp$ma_bmd[, 1]), 1],
      method = "monoH.FC", ties=mean
    )
    temp$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
    t_names <- names(temp)

    idx <- grep("Fitted_Model", t_names)
    for (ii in idx) {
      # temp[[ii]]$prior <- priors[[which(ii == idx)]]
      temp[[ii]]$data <- data
      temp[[ii]]$model <- tolower(trimws(gsub("Model: ", "", temp[[ii]]$full_model)))
      if (temp[[ii]]$model == "quantal-linear") {
        temp[[ii]]$model <- "qlinear"
      }

      data_temp <- temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[, 1]), ]
      if (length(data_temp) > 0) {
        data_temp <- data_temp[!is.na(data_temp[, 1]), ]
        if (nrow(data_temp) > 6) {
          te <- splinefun(sort(data_temp[, 2, drop = F]), sort(data_temp[, 1, drop = F]), method = "monoH.FC",ties=mean)
          temp[[ii]]$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
          if (max(data_temp[, 2]) < 1 - alpha) {
            temp[[ii]]$bmd[3] <- 1e300
          }
        } else {
          temp[[ii]]$bmd <- c(NA, NA, NA)
        }
      }else{
        temp[[ii]]$bmd <- c(NA,NA,NA)
      }

      #te <- splinefun(temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[, 1]), 2], temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[, 1]), 1], method = "monoH.FC",ties=mean)
      #temp[[ii]]$bmd <- c(te(0.5), te(alpha), te(1 - alpha))

      #add NAs if bad hessian or NaN BMD
      if(det(temp[[ii]]$covariance) < 0 || is.na(temp[[ii]]$bmd[1])){
        temp$posterior_probs[ii] <- NA
      }

      names(temp[[ii]]$bmd) <- c("BMD", "BMDL", "BMDU")
      # names(temp)[ii] <- sprintf("Individual_Model_%s", ii)
      names(temp)[ii] <- sprintf("Indiv_%s_%s", model_list[ii], distribution_list[ii])

      tmp_id <- which(names(temp) == "BMD_CDF")
      #  temp = temp[-tmp_id]
      temp[[ii]]$options <- c(o1, o2)

      #names(temp$posterior_probs) <- paste(model_list, distribution_list, sep="_")

    }
    names(temp$posterior_probs) <- paste(model_list, distribution_list, sep="_")

    class(temp) <- c("BMDdichotomous_MA", "BMDdichotomous_MA_laplace")
  } else {
    # MCMC run
    temp_r <- .run_ma_dichotomous(
      data, priors, model_i,
      model_p, TRUE, o1, o2
    )
    tempn <- temp_r$ma_results
    tempm <- temp_r$mcmc_runs
    # clean up the run
    idx <- grep("Fitted_Model", names(tempn))
    temp <- list()
    jj <- 1
    for (ii in idx) {
      temp[[jj]] <- list()
      temp[[jj]]$mcmc_result <- tempm[[ii]]
      temp[[jj]]$full_model <- tempn[[ii]]$full_model
      temp[[jj]]$parameters <- tempn[[ii]]$parameters
      temp[[jj]]$covariance <- tempn[[ii]]$covariance
      temp[[jj]]$maximum <- tempn[[ii]]$maximum
      temp[[jj]]$bmd_dist <- tempn[[ii]]$bmd_dist

      # temp[[jj]]$prior <- priors[[which(ii == idx)]]
      temp[[jj]]$data <- data
      temp[[jj]]$model <- tolower(trimws(gsub("Model: ", "", tempn[[ii]]$full_model)))
      temp[[jj]]$options <- c(o1, o2)
      if (temp[[jj]]$model == "quantal-linear") {
        temp[[jj]]$model <- "qlinear"
      }

      data_temp <- temp[[jj]]$bmd_dist[!is.infinite(temp[[jj]]$bmd_dist[, 1]), ]
      if (length(data_temp) > 0) {
        data_temp <- data_temp[!is.na(data_temp[, 1]), ]
        if (nrow(data_temp) > 6) {
          te <- splinefun(sort(data_temp[, 2, drop = F]), sort(data_temp[, 1, drop = F]), method = "monoH.FC",ties=mean)
          temp[[ii]]$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
          if (max(data_temp[, 2]) < 1 - alpha) {
            temp[[ii]]$bmd[3] <- 1e300
          }
        } else {
          temp[[ii]]$bmd <- c(NA, NA, NA)
        }
      }else{
        temp[[ii]]$bmd <- c(NA,NA,NA)
      }

      #te <- splinefun(temp[[jj]]$bmd_dist[!is.infinite(temp[[jj]]$bmd_dist[, 1]), 2], temp[[jj]]$bmd_dist[!is.infinite(temp[[jj]]$bmd_dist[, 1]), 1], method = "monoH.FC",ties=mean)
      #temp[[jj]]$bmd <- c(te(0.5), te(alpha), te(1 - alpha))

      #add NAs if bad hessian or NaN BMD
      if(det(temp[[ii]]$covariance) < 0 || is.na(temp[[ii]]$bmd[1])){
        tempn$posterior_probs[ii] <- NA
      }

      class(temp[[jj]]) <- "BMDdich_fit_MCMC"
      jj <- jj + 1
    }
    # for (ii in idx_mcmc)
    # names(temp) <- sprintf("Individual_Model_%s", 1:length(priors))
    names(temp) <- sprintf("Indiv_%s_%s", model_list, distribution_list)

    temp$ma_bmd <- tempn$BMD_CDF
    te <- splinefun(temp$ma_bmd[!is.infinite(temp$ma_bmd[, 1]), 2], temp$ma_bmd[!is.infinite(temp$ma_bmd[, 1]), 1], method = "monoH.FC",ties=mean)
    temp$bmd <- c(te(0.5), te(alpha), te(1 - alpha))
    temp$posterior_probs <- tempn$posterior_probs
    names(temp$posterior_probs) <- paste(model_list, distribution_list, sep="_")
    #temp$post_prob
    class(temp) <- c("BMDdichotomous_MA", "BMDdichotomous_MA_mcmc")
  }
  names(temp$bmd) <- c("BMD", "BMDL", "BMDU")
  return(temp)
}
