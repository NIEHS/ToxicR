# Copyright 2022  NIEHS <matt.wheeler@nih.gov>
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


.evaluate_alpha <- function(...) {
  args <- list(...)

  cl <- match.call()
  ev <- match("alpha", names(cl))
  if (is.na(ev)) {
    alpha <- 0.05
  } else {
    alpha <- args$alpha
    if ((alpha <= 0) || alpha > 0.5) {
      stop("alpha must be in (0,0.5]")
    }
  }

  return(alpha)
}


.summary_continuous_max <- function(object, ...) {
  model <- object
  returnV <- list()
  #alpha <- .evaluate_alpha(...)
  alpha <- model$options[4]

  if (is.null(model$prior)) {
    returnV$fit_method <- "MLE"
    returnV$prior <- NA
  } else {
    returnV$fit_method <- "Bayesian:MAP"
    returnV$prior <- model$prior
  }
  returnV$fit <- model$full_model

  temp_function <- tryCatch(splinefun(model$bmd_dist[,2],model$bmd_dist[,1],method="monoH.FC"),error = function(e) {NULL})
  if (is.null(temp_function)){
     returnV$BMD <- c(NA,NA,NA)
  }else{
      returnV$BMD <- temp_function(1-c(1-alpha,0.5,alpha))
  }

  names(returnV$BMD) <- c("BMDL", "BMD", "BMDU")
  returnV$alpha <- alpha

  temp <- model$Deviance - matrix(rep(model$Deviance[5, ], 5), nrow = 5, ncol = 2, byrow = T)
  if (grepl("NCV", model$full_model)) {
    temp <- rbind(temp[3, ], temp[2, ])
  } else {
    temp <- rbind(temp[1, ], temp[2, ])
  }
  temp[, 1] <- -1 * temp[, 1]
  if(all(is.infinite(temp[,1]))){
    stop("Summary failed. Need multiple observations per dose.")
  }
  returnV$GOF <- cbind(temp, pchisq(temp[, 1], temp[, 2], lower.tail = F))
  colnames(returnV$GOF) <- c("-2LL", "DF", "P-Value")
  class(returnV) <- "summary_continuous_max"
  return(returnV)
}

.print_summary_continuous_max <- function(x, ...) { # nolint
  s_fit <- x

  if (grepl("MLE", s_fit$fit_method)) {
    cat(sprintf("Summary of single model fit (%s) using ToxicR\n", "MLE"))
    cat(s_fit$fit, "\n")
  } else {
    cat(sprintf("Summary of single model fit (%s) using ToxicR\n\n", "Bayesian-MAP"))
    s_fit$GOF[, 2] <- round(s_fit$GOF[, 2], 2)
  }
  cat("\n")

  cat("BMD: ")
  cat(sprintf("%1.2f (%1.2f, %1.2f) %1.1f%% CI\n", s_fit$BMD[2], s_fit$BMD[1], s_fit$BMD[3], 100 * (1 - 2 * s_fit$alpha)))
  cat("\n")
  cat("Model GOF\n")
  cat("--------------------------------------------------\n")
  s_fit$GOF[, 1] <- round(s_fit$GOF[, 1], 2)

  s_fit$GOF[, 3] <- round(s_fit$GOF[, 3], 3)
  rownames(s_fit$GOF) <- c("Test: Mean Adequate", "Test: Mean/Variance Adequate")
  print(s_fit$GOF)
}


.summary_continuous_mcmc <- function(object, ...) {
  model <- object
  returnV <- list()

  #alpha <- .evaluate_alpha(...)
  #alpha is different for dichotomous
  if(any(grepl("dich", class(model)))){
    alpha <- model$options[2]
  }else{
    alpha <- model$options[4]
  }

  returnV$fit_method <- "Bayesian:MCMC"
  returnV$prior <- model$prior
  returnV$fit <- model$full_model

  temp_function <- tryCatch(splinefun(model$bmd_dist[,2],model$bmd_dist[,1],method="monoH.FC"),error = function(e) {NULL})
  if (is.null(temp_function)){
     returnV$BMD <- c(NA,NA,NA)
  }else{
      returnV$BMD <- temp_function(1-c(1-alpha,0.5,alpha))
  }
  names(returnV$BMD) <- c("BMDL", "BMD", "BMDU")
  returnV$alpha <- alpha
  returnV$eff_size <- coda::effectiveSize(model$mcmc_result$BMD_samples)
  returnV$geweke_z <- coda::geweke.diag(coda::as.mcmc(model$mcmc_result$BMD_samples), frac1 = 0.3, frac2 = 0.4)$z
  class(returnV) <- "summary_mcmc"
  return(returnV)
}

.print_summary_continuous_mcmc <- function(x, ...) {
  s_fit <- x

  cat(sprintf("Summary of single model fit (%s) using ToxicR\n", "MCMC"))
  cat(s_fit$fit, "\n")

  cat("\n")

  cat("BMD: ")
  cat(sprintf("%1.2f (%1.2f, %1.2f) %1.1f%% CI\n", s_fit$BMD[2], s_fit$BMD[1], s_fit$BMD[3], 100 * (1 - 2 * s_fit$alpha)))
  cat("\n")
  cat("Convergence Diagnostics on BMD\n")
  cat("--------------------------------------------------\n")
  cat(sprintf("Effective Sample Size: %1.2f\n\n", s_fit$eff_size))
  cat(sprintf(
    "Geweke Z-score that mean of first 30%% of \nMCMC chain is different from last 40%%\nZ-Score: %1.3f  P-value %1.3f\n",
    s_fit$geweke_z, 2 * pnorm(abs(s_fit$geweke_z), lower.tail = F)
  ))
}

.summary_ma_max <- function(object, ...) {
  model <- object
  #alpha <- .evaluate_alpha(...)
  tmp_idx <- grep("Indiv_", names(model))
  #alpha is different for dichotomous
  if(any(grepl("dichotomous", class(model)))){
    alpha <- model[[tmp_idx[1]]]$options[2]
  }else{
    alpha <- model[[tmp_idx[1]]]$options[4]
  }

  returnV <- list()

  returnV$fit_method <- "Bayesian:MAX"
  returnV$fit_table <- data.frame(post_p = round(model$posterior_probs, 3))

  temp_mfit <- rep(" ", length(tmp_idx)) # model name
  temp_BMD <- rep(" ", length(tmp_idx)) # bmd
  temp_BMDL <- rep(" ", length(tmp_idx)) # bmdl
  temp_BMDU <- rep(" ", length(tmp_idx)) # bmdu

  for (ii in tmp_idx) {
    tmp_fit <- model[[ii]]
    data_temp <- tmp_fit$bmd_dist
    dist <- data_temp[!is.infinite(data_temp[, 1]) & !is.na(data_temp[, 1]), ]
    dist <- data_temp[!is.nan(data_temp[, 1])]
    if (length(dist) > 10 & !identical(dist, numeric(0))) {
      temp_function <- splinefun(data_temp[, 2], data_temp[, 1], method = "monoH.FC",ties=mean)
      temp_bmds <- temp_function(1 - c(1 - alpha, 0.5, alpha))
      temp_mfit[ii] <- sub("Model: ", "", tmp_fit$full_model)
      temp_BMD[ii] <- round(temp_bmds[2], 3)
      temp_BMDL[ii] <- round(temp_bmds[1], 3)
      temp_BMDU[ii] <- round(temp_bmds[3], 3)
    } else {
      temp_mfit[ii] <- sub("Model: ", "", tmp_fit$full_model)
      temp_BMD[ii] <- NA
      temp_BMDL[ii] <- NA
      temp_BMDU[ii] <- NA
    }
  }

  returnV$fit_table$model_names <- temp_mfit
  returnV$fit_table$BMD <- temp_BMD
  returnV$fit_table$BMDL <- temp_BMDL
  returnV$fit_table$BMDU <- temp_BMDU

  tmp_idx <- order(returnV$fit_table$post_p, decreasing = T)
  returnV$fit_table <- returnV$fit_table[tmp_idx, c(2, 3, 4, 5, 1)]

  warnFunc <- function(w) {
   return(-2)
  }
  errorFunc <- function(w){
    return(-2)
  }
  ma_temp_function = tryCatch(
     {
         splinefun(model$ma_bmd[, 2], model$ma_bmd[, 1], method = "monoH.FC",ties=mean)
     },
      warning = warnFunc
      ,
      error = errorFunc
  )
  if (is.numeric(ma_temp_function)){
    returnV$BMD <- c(NA,NA,NA)
  }
  
  returnV$BMD <- ma_temp_function(1 - c(1 - alpha, 0.5, alpha))
  names(returnV$BMD) <- c("BMDL", "BMD", "BMDU")
  returnV$alpha <- alpha

  class(returnV) <- "ma_summary_max"
  return(returnV)
}


.summary_ma_mcmc <- function(object, ...) {
  model <- object
  #alpha <- .evaluate_alpha(...)
  tmp_idx <- grep("Indiv_", names(model))
  #alpha is different for dichotomous
  if(any(grepl("dichotomous", class(model)))){
    alpha <- model[[tmp_idx[1]]]$options[2]
  }else{
    alpha <- model[[tmp_idx[1]]]$options[4]
  }

  returnV <- list()

  returnV$fit_method <- "Bayesian:MCMC"
  returnV$fit_table <- data.frame(post_p = round(model$posterior_probs, 3))

  temp_mfit <- rep(" ", length(tmp_idx)) # model name
  temp_BMD <- rep(" ", length(tmp_idx)) # bmd
  temp_BMDL <- rep(" ", length(tmp_idx)) # bmdl
  temp_BMDU <- rep(" ", length(tmp_idx)) # bmdu

  for (ii in tmp_idx) {
    tmp_fit <- model[[ii]]
    data_temp <- tmp_fit$bmd_dist
    dist <- data_temp[!is.infinite(data_temp[, 1]) & !is.na(data_temp[, 1]), ]
    dist <- data_temp[!is.nan(data_temp[, 1])]
    if (length(dist) > 10 & !identical(dist, numeric(0))) {
      temp_function <- splinefun(data_temp[, 2], data_temp[, 1], method = "monoH.FC",ties=mean)
      temp_bmds <- temp_function(1 - c(1 - alpha, 0.5, alpha))
      temp_mfit[ii] <- sub("Model: ", "", tmp_fit$full_model)
      temp_BMD[ii] <- round(temp_bmds[2], 3)
      temp_BMDL[ii] <- round(temp_bmds[1], 3)
      temp_BMDU[ii] <- round(temp_bmds[3], 3)
    } else {
      temp_mfit[ii] <- sub("Model: ", "", tmp_fit$full_model)
      temp_BMD[ii] <- NA
      temp_BMDL[ii] <- NA
      temp_BMDU[ii] <- NA
    }
  }

  returnV$fit_table$model_names <- temp_mfit
  returnV$fit_table$BMD <- temp_BMD
  returnV$fit_table$BMDL <- temp_BMDL
  returnV$fit_table$BMDU <- temp_BMDU

  tmp_idx <- order(returnV$fit_table$post_p, decreasing = T)
  returnV$fit_table <- returnV$fit_table[tmp_idx, c(2, 3, 4, 5, 1)]

  warnFunc <- function(w) {
    return(-2)
  }
  errorFunc <- function(w){
    return(-2)
  }
  ma_temp_function = tryCatch(
    {
      splinefun(model$ma_bmd[, 2], model$ma_bmd[, 1], method = "monoH.FC",ties=mean)
    },
    warning = warnFunc
    ,
    error = errorFunc
  )
  if (is.numeric(ma_temp_function)){
    returnV$BMD <- c(NA,NA,NA)
  }
  

  returnV$BMD <- ma_temp_function(1 - c(1 - alpha, 0.5, alpha))
  names(returnV$BMD) <- c("BMDL", "BMD", "BMDU")
  returnV$alpha <- alpha

  class(returnV) <- "ma_summary_mcmc"
  return(returnV)
}

.print_summary_ma <- function(x, ...) { # nolint
  s_fit <- x
  cat("Summary of single MA BMD\n\n")
  cat("Individual Model BMDS\n")
  cat(paste("Model", strrep(" ", 59), sep = ""), "\t\t BMD (BMDL, BMDU)\tPr(M|Data)\n")
  cat("___________________________________________________________________________________________\n")
  badd <- c()
  for (ii in 1:nrow(s_fit$fit_table)) {
    tmp_length <- nchar(s_fit$fit_table[ii, 1])
    # pad <- paste(substr(s_fit$fit_table[ii, 1], 1, 38), strrep(" ", 39 - tmp_length), sep = "")
    pad <- paste(substr(s_fit$fit_table[ii, 1], 1, 61), strrep(" ", 62 - tmp_length), sep = "")
    if(all(is.na(s_fit$fit_table[ii,2:5]))){
      badd <- c(badd, ii)
    }else{
      cat(sprintf(
        "%s\t\t\t%1.2f (%1.2f ,%1.2f) \t %1.3f\n", pad, as.numeric(s_fit$fit_table[ii, 2]),
        as.numeric(s_fit$fit_table[ii, 3]), as.numeric(s_fit$fit_table[ii, 4]), as.numeric(s_fit$fit_table[ii, 5])
      ))
    }
  }
  #print at end the values with NA BMD BMDL BMDU
  for(jj in badd){
    tmp_length <- nchar(s_fit$fit_table[jj, 1])
    pad <- paste(substr(s_fit$fit_table[jj, 1], 1, 61), strrep(" ", 62 - tmp_length), sep = "")
    cat(sprintf(
      "%s\t\t\t%1.2f (%1.2f ,%1.2f) \t \t %1.3f\n", pad, as.numeric(s_fit$fit_table[jj, 2]),
      as.numeric(s_fit$fit_table[jj, 3]), as.numeric(s_fit$fit_table[jj, 4]), as.numeric(s_fit$fit_table[jj, 5])
    ))
  }
  cat("___________________________________________________________________________________________\n")
  
  cat("Model Average BMD: ")
  cat(sprintf(
    "%1.2f (%1.2f, %1.2f) %1.1f%% CI\n", s_fit$BMD[2], s_fit$BMD[1], s_fit$BMD[3],
    100 * (1 - 2 * s_fit$alpha)
  ))
}
