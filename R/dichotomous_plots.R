# Dichotomous functions are defined here
{
  .logit <- function(p) {
    return(log(p / (1 - p)))
  }

  # dichotomous hill
  .dich_hill_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    n <- 1 / (1 + exp(-parms[2]))
    a <- parms[3]
    b <- parms[4]
    rval <- g + (1 - g) * n * (1 / (1 + exp(-a - b * log(d))))
    return(rval)
  }
  # dichotomous log-logistic
  .dich_llogist_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    a <- parms[2]
    b <- parms[3]
    rval <- g + (1 - g) * (1 / (1 + exp(-a - b * log(d))))
    return(rval)
  }
  # dichotomous log-probit
  .dich_lprobit_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    a <- parms[2]
    b <- parms[3]
    rval <- g + (1 - g) * (1-pnorm(-a - b * log(d)))
    return(rval)
  }

  # dichotomous weibull
  .dich_weibull_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    a <- parms[2]
    b <- parms[3]
    rval <- g + (1 - g) * (1 - exp(-b * d^a))
    return(rval)
  }

  # dichotomous gamma
  .dich_gamma_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    a <- parms[2]
    b <- parms[3]
    rval <- g + (1 - g) * pgamma(b * d, a, 1)
    return(rval)
  }

  # dichtomous logistic
  .dich_logist_f <- function(parms, d) {
    rval <- 1 / (1 + exp(-parms[1] - parms[2] * d))
    return(rval)
  }

  # dichtomous probit
  .dich_probit_f <- function(parms, d) {
    rval <- pnorm(parms[1] + parms[2] * d)
    return(rval)
  }

  .dich_qlinear_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    a <- parms[2]
    return(g + (1 - g) * (1 - exp(-a * d)))
  }

  .dich_multistage_f <- function(parms, d) {
    g <- 1 / (1 + exp(-parms[1]))
    rval <- d * 0
    for (ii in 2:length(parms)) {
      rval <- rval - parms[ii] * d^(ii - 1)
    }
    return(g + (1 - g) * (1 - exp(rval)))
  }
}

{
  .plot.BMDdich_fit_MCMC <- function(x, ...) {
    fit <- x
    temp_args <- list(...)

    if (!exists("qprob", temp_args)) {
      qprob <- 0.05
    } else {
      qprob <- temp_args$qprob
    }
    density_col <- "red"
    credint_col <- "azure2"
    BMD_DENSITY <- T

    if (qprob < 0 || qprob > 0.5) {
      stop("Quantile probability must be between 0 and 0.5")
    }

    # How this is calculated?
    # This part - how it is derived?
    probs <- (0.5 + fit$data[, 2, drop = T]) / (1.0 + fit$data[, 3, drop = T])
    se <- sqrt(probs * (1 - probs) / fit$data[, 3, drop = T])


    doses <- fit$data[, 1, drop = T]
    uerror <- apply(cbind(probs * 0 + 1, probs + se), 1, min, na.rm = TRUE)
    lerror <- apply(cbind(probs * 0, probs - se), 1, max, na.rm = TRUE)

    dose <- c(doses, doses)
    Response <- c(uerror, lerror)

    # Basic structure of display
    # main should show the models' information, I think this part should be fixed.

    # Dichotomous's response is between 0 to 1
    # Change this to ggplot object


    # plot(dose,Response,type='n',main=fit$full_model)
    # We need to adjust the range here too
    # S3 object not fitted here for the title part

    test_doses <- seq(min(doses), max(doses) * 1.03, (max(doses) * 1.03 - min(doses)) / 100)

    if (fit$model == "hill") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_hill_f, d = test_doses)
    }
    if (fit$model == "gamma") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_gamma_f, d = test_doses)
    }
    if (fit$model == "logistic") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_logist_f, d = test_doses)
    }
    if (fit$model == "log-logistic") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_llogist_f, d = test_doses)
    }
    if (fit$model == "probit") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_probit_f, d = test_doses)
    }
    if (fit$model == "log-probit") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_lprobit_f, d = test_doses)
    }
    if (fit$model == "multistage") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_multistage_f, d = test_doses)
    }
    if (fit$model == "qlinear") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_qlinear_f, d = test_doses)
    }

    if (fit$model == "weibull") {
      Q <- apply(fit$mcmc_result$PARM_samples, 1, .dich_weibull_f, d = test_doses)
    }


    temp <- fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp <- temp[!is.infinite(temp)]
    test <- density(temp)


    Q <- t(Q)

    me <- colMeans(Q, , na.rm = TRUE)
    lq <- apply(Q, 2, quantile, probs = qprob, na.rm = TRUE)
    uq <- apply(Q, 2, quantile, probs = 1 - qprob, na.rm = TRUE)

    plot_gg <- ggplot() +
      geom_errorbar(aes(x = doses, ymin = lerror, ymax = uerror), color = "grey") +
      labs(x = "Dose", y = "Proportion", title = paste(fit$full_model, sep = ",  Fit Type: ")) +
      theme_minimal() +
      xlim(0 - 5 * max(test_doses), 5 * max(test_doses))

    # Spline function is used to test column average from MCMC

    temp_fit <- splinefun(test_doses, me)

    # Object 2
    # Polygon changed to Geom_ribbon

    plot_gg <- plot_gg + geom_ribbon(aes(x = test_doses, ymin = lq, ymax = uq), fill = "blue", alpha = 0.1)
    plot_gg <- plot_gg +
      geom_line(aes(x = test_doses, y = me), col = "blue", linewidth = 2) + geom_point(aes(x = doses, y = probs))


    plot_gg <- plot_gg +
      geom_segment(aes(
        x = fit$bmd[2], y = temp_fit(fit$bmd[1]), xend = fit$bmd[3],
        yend = temp_fit(fit$bmd[1])
      ), color = "darkslategrey", linewidth = 1.2, alpha = 0.9) +
      annotate(
        geom = "text", x = fit$bmd[2], y = temp_fit(fit$bmd[1]),
        label = "[", size = 10, color = "darkslategrey", alpha = 0.9
      ) +
      annotate(
        geom = "text", x = fit$bmd[3], y = temp_fit(fit$bmd[1]),
        label = "]", size = 10, color = "darkslategrey", alpha = 0.9
      ) +
      annotate(
        geom = "point", x = fit$bmd[1], y = temp_fit(fit$bmd[1]),
        size = 5, color = "darkslategrey", shape = 17, alpha = 0.9
      )


    # plot_gg<-plot_gg+
    #   geom_segment(aes(x=fit$bmd, y=temp_fit(x=fit$bmd), xend=fit$bmd, yend=min(Response,me)*0.95),color="Red")
    #


    # Adding density

    # Object 3 - Density object - In Shiny we can on / off this
    # Density - Polygon/Other option?
    if (BMD_DENSITY == TRUE) {
      Dens <- density(temp, cut = c(5 * max(test_doses)), n = 1000, from = 0, to = max(test_doses), na.rm = TRUE)

      Dens$y <- Dens$y / max(Dens$y) * (max(uerror) - min(lerror)) * 0.6
      temp <- which(Dens$x < max(test_doses))
      D1_y <- Dens$y[temp]
      D1_x <- Dens$x[temp]
      qm <- min(lerror)
      scale <- (max(uerror) - min(lerror)) / max(D1_y) * .40


      plot_gg <- plot_gg +
        geom_polygon(aes(
          x = c(max(0, min(D1_x)), D1_x, max(D1_x)),
          y = c(min(lerror), min(lerror) + D1_y * scale, min(lerror))
        ),
        fill = "blueviolet", alpha = 0.6
        )
    }

    return(plot_gg + coord_cartesian(xlim = c(min(doses), max(doses)), expand = F))
  }

  .plot.BMDdich_fit_maximized <- function(x, ...) {
    fit <- x
    temp_args <- list(...)
    if (!exists("qprob", temp_args)) {
      qprob <- 0.05
    } else {
      qprob <- temp_args$qprob
    }


    density_col <- "red"
    credint_col <- "azure2"

    probs <- (0.5 + fit$data[, 2, drop = T]) / (1.0 + fit$data[, 3, drop = T])
    se <- sqrt(probs * (1 - probs) / fit$data[, 3, drop = T])


    doses <- fit$data[, 1, drop = T]
    uerror <- apply(cbind(probs * 0 + 1, probs + se), 1, min)
    lerror <- apply(cbind(probs * 0, probs - se), 1, max)

    dose <- c(doses, doses)
    Response <- c(uerror, lerror)



    test_doses <- seq(min(doses), max(doses) * 1.03, (max(doses) * 1.03 - min(doses)) / 100)


    # Need to check loop
    if (fit$model == "hill") {
      me <- .dich_hill_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "gamma") {
      me <- .dich_gamma_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "logistic") {
      me <- .dich_logist_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "log-logistic") {
      me <- .dich_llogist_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "probit") {
      me <- .dich_probit_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "log-probit") {
      me <- .dich_lprobit_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "multistage") {
      me <- .dich_multistage_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "qlinear") {
      me <- .dich_qlinear_f(fit$parameters, d = test_doses)
    }
    if (fit$model == "weibull") {
      me <- .dich_weibull_f(fit$parameters, d = test_doses)
    }



    temp_fit <- splinefun(test_doses, me)


    plot_gg <- ggplot() +
      geom_errorbar(aes(x = doses, ymin = lerror, ymax = uerror), color = "grey") +
      xlim(c(min(dose) - 5 * max(dose), max(dose) * 5)) +
      labs(x = "Dose", y = "Proportion", title = paste(fit$full_model, sep = ",  Fit Type: ")) +
      theme_minimal()



    # BMD Estimates fit - MLE/Laplace why they don't have it yet..?


    plot_gg <- plot_gg +
      geom_line(aes(x = test_doses, y = me), col = "blue", linewidth = 1.2) + geom_point(aes(x = doses, y = probs))


    plot_gg <- plot_gg +
      geom_segment(aes(
        x = fit$bmd[2], y = temp_fit(fit$bmd[1]), xend = fit$bmd[3],
        yend = temp_fit(fit$bmd[1])
      ), color = "darkslategrey", linewidth = 1.2, alpha = 0.9) +
      annotate(
        geom = "text", x = fit$bmd[2], y = temp_fit(fit$bmd[1]),
        label = "[", size = 10, color = "darkslategrey", alpha = 0.9
      ) +
      annotate(
        geom = "text", x = fit$bmd[3], y = temp_fit(fit$bmd[1]),
        label = "]", size = 10, color = "darkslategrey", alpha = 0.9
      ) +
      annotate(
        geom = "point", x = fit$bmd[1], y = temp_fit(fit$bmd[1]),
        size = 5, color = "darkslategrey", shape = 17, alpha = 0.9
      )

    return(plot_gg + coord_cartesian(xlim = c(min(doses), max(doses)), expand = F))
  }

  .plot.BMDdichotomous_MA <- function(x, ...) {
    A <- x
    model_no <- x_axis <- y_axis <- cols <- NULL
    temp_args <- list(...)

    if (!exists("qprob", temp_args)) {
      qprob <- 0.05
    } else {
      qprob <- temp_args$qprob
    }

    density_col <- "blueviolet"
    credint_col <- "azure2"
    fit_origin <- A # Updated SL
    class_list <- names(A)
    fit_idx <- grep("Indiv_", class_list) # 06/18/21 SL


    # plot the model average curve
    if ("BMDdichotomous_MA_mcmc" %in% class(A)) { # mcmc run
      A$posterior_probs[!is.finite(A$posterior_probs)] = 0
      n_samps <- nrow(A[[fit_idx[1]]]$mcmc_result$PARM_samples)
      data_d <- A[[fit_idx[1]]]$data
      max_dose <- max(data_d[, 1])
      min_dose <- min(data_d[, 1])
      test_doses <- seq(min_dose, max_dose, (max_dose - min_dose) / 500)
      ma_samps <- sample(fit_idx, n_samps, replace = TRUE, prob = unname(A$posterior_probs))
      temp_f <- matrix(0, n_samps, length(test_doses))
      temp_bmd <- rep(0, length(test_doses))



      probs <- (0.5 + data_d[, 2, drop = T]) / (1.0 + data_d[, 3, drop = T])
      se <- sqrt(probs * (1 - probs) / data_d[, 3, drop = T])
      doses <- data_d[, 1, drop = T]
      uerror <- apply(cbind(probs * 0 + 1, probs + se), 1, min)
      lerror <- apply(cbind(probs * 0, probs - se), 1, max)

      dose <- c(doses, doses)
      Response <- c(uerror, lerror)

      plot_gg <- ggplot() +
        geom_errorbar(aes(x = doses, ymin = lerror, ymax = uerror), color = "grey") +
        xlim(c(-5 * max(dose)), 5 * max(dose)) +
        labs(x = "Dose", y = "Proportion", title = "Model : Dichotomous MA") +
        theme_minimal()


      for (ii in 1:n_samps) {
        fit <- A[[fit_idx[ma_samps[ii]]]]

        if (fit$model == "hill") {
          temp_f[ii, ] <- .dich_hill_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "gamma") {
          temp_f[ii, ] <- .dich_gamma_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "logistic") {
          temp_f[ii, ] <- .dich_logist_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "log-logistic") {
          temp_f[ii, ] <- .dich_llogist_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "probit") {
          temp_f[ii, ] <- .dich_probit_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "log-probit") {
          temp_f[ii, ] <- .dich_lprobit_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "multistage") {
          temp_f[ii, ] <- .dich_multistage_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "qlinear") {
          temp_f[ii, ] <- .dich_qlinear_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
        if (fit$model == "weibull") {
          temp_f[ii, ] <- .dich_weibull_f(fit$mcmc_result$PARM_samples[ii, ], test_doses)
          temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
        }
      }


      me <- colMeans(temp_f) # Why col means instead of median? check line 372 for continues_plots.R
      # me <- apply(temp_f,2,quantile, probs = 0.5,na.rm = TRUE) # BMD
      lq <- apply(temp_f, 2, quantile, probs = qprob, na.rm = TRUE)
      uq <- apply(temp_f, 2, quantile, probs = 1 - qprob, na.rm = TRUE)

      plot_gg <- plot_gg +
        geom_ribbon(aes(x = test_doses, ymin = lq, ymax = uq), fill = "blue", alpha = 0.1)

      plot_gg <- plot_gg +
        geom_line(aes(x = test_doses, y = me), col = "blue", linewidth = 2) +
        geom_point(aes(x = doses, y = probs))


      temp_fit <- splinefun(test_doses, me)

      fit <- A

      plot_gg <- plot_gg +
        geom_segment(aes(
          x = A$bmd[2], y = temp_fit(A$bmd[1]), xend = A$bmd[3],
          yend = temp_fit(A$bmd[1])
        ), color = "darkslategrey", linewidth = 1.2, alpha = 0.9)
      plot_gg <- plot_gg +
        annotate(
          geom = "text", x = A$bmd[2], y = temp_fit(A$bmd[1]),
          label = "[", size = 10, color = "darkslategrey", alpha = 0.9
        ) +
        annotate(
          geom = "text", x = A$bmd[3], y = temp_fit(A$bmd[1]),
          label = "]", size = 10, color = "darkslategrey", alpha = 0.9
        ) +
        annotate(
          geom = "point", x = A$bmd[1], y = temp_fit(A$bmd[1]),
          size = 5, color = "darkslategrey", shape = 17, alpha = 0.9
        )


      # Density needs to be re derived ... based on the continous logic in the MA case
      temp <- temp_bmd[!is.nan(temp_bmd)]
      temp <- temp[!is.infinite(temp)]
      temp <- temp[temp < 5 * max(doses)]

      Dens <- density(temp, cut = c(5 * max(test_doses)), n = 1000, from = 0, to = max(test_doses))

      Dens$y <- Dens$y / max(Dens$y) * (max(uerror) - min(lerror)) * 0.6
      temp <- which(Dens$x < max(test_doses))
      D1_y <- Dens$y[temp]
      D1_x <- Dens$x[temp]
      qm <- min(lerror)
      scale <- (max(uerror) - min(lerror)) / max(D1_y) * .40


      plot_gg <- plot_gg +
        geom_polygon(aes(
          x = c(max(0, min(D1_x)), D1_x, max(D1_x)),
          y = c(min(lerror), min(lerror) + D1_y * scale, min(lerror))
        ),
        fill = "blueviolet", alpha = 0.6
        )
      # geom_polygon(aes(x=c(0,D1_x,max(doses)),y=c(qm,qm+D1_y,qm)), fill = "blueviolet", alpha=0.6)


      # plot the individual models proportional to their weight

      # Reset plot cage
      # temp_f <- rep(0,length(test_doses))
      temp_house <- matrix(nrow = length(fit_idx), ncol = length(temp_f))


      df <- NULL

      for (ii in 1:length(fit_idx)) {
        if (!is.finite(A$posterior_probs[ii])){
           A$posterior_probs[ii] = 0
        }
        if (A$posterior_probs[ii] > 0.05) {
          fit <- A[[fit_idx[ii]]]
          if (fit$model == "hill") {
            f <- .dich_hill_f(fit$parameters, test_doses)
          }
          if (fit$model == "gamma") {
            f <- .dich_gamma_f(fit$parameters, test_doses)
          }
          if (fit$model == "logistic") {
            f <- .dich_logist_f(fit$parameters, test_doses)
          }
          if (fit$model == "log-logistic") {
            f <- .dich_llogist_f(fit$parameters, test_doses)
          }
          if (fit$model == "probit") {
            f <- .dich_probit_f(fit$parameters, test_doses)
          }
          if (fit$model == "log-probit") {
            f <- .dich_lprobit_f(fit$parameters, test_doses)
          }
          if (fit$model == "multistage") {
            f <- .dich_multistage_f(fit$parameters, test_doses)
          }
          if (fit$model == "qlinear") {
            f <- .dich_qlinear_f(fit$parameters, test_doses)
          }
          if (fit$model == "weibull") {
            f <- .dich_weibull_f(fit$parameters, test_doses)
          }

          col <- "coral3"
          temp_df <- data.frame(x_axis = test_doses, y_axis = f, cols = col, model_no = ii, alpha_lev = unname(A$posterior_probs[ii]))
          df <- rbind(df, temp_df)

          # SL Updated 06/18/21 -- Transparency update based on posterior probability and Y scale for dichotomous case
          temp_data <- df %>%
            filter(model_no == ii)

          plot_gg <- plot_gg +
            geom_line(data = temp_data, aes(x = x_axis, y = y_axis, color = cols), alpha = unique(temp_data$alpha_lev), show.legend = F) +
            theme_minimal()
        }
      }


      return(plot_gg + coord_cartesian(xlim = c(min(doses), max(doses)), expand = F))
    } else if ("BMDdichotomous_MA_laplace" %in% class(A)) { # mcmc run
      A$posterior_probs[!is.finite(A$posterior_probs)] = 0
      class_list <- names(A)
      fit_idx <- grep("Indiv_", class_list)
      num_model <- length(A$posterior_probs)

      data_d <- A[[1]]$data
      max_dose <- max(data_d[, 1])
      min_dose <- min(data_d[, 1])
      test_doses <- seq(min_dose, max_dose, (max_dose - min_dose) / 500)

      # Create 0 matrix
      temp_f <- matrix(0, num_model, length(test_doses))
      temp_bmd <- rep(0, length(test_doses))

      probs <- (0.5 + data_d[, 2, drop = T]) / (1.0 + data_d[, 3, drop = T])
      se <- sqrt(probs * (1 - probs) / data_d[, 3, drop = T])
      doses <- data_d[, 1, drop = T]
      uerror <- apply(cbind(probs * 0 + 1, probs + se), 1, min)
      lerror <- apply(cbind(probs * 0, probs - se), 1, max)

      dose <- c(doses, doses)
      Response <- c(uerror, lerror)

      # Line plot for based on each cases
      for (ii in 1:num_model) {
        fit_loop <- A[[ii]]

        if (fit_loop$model == "hill") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_hill_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "gamma") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_gamma_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "logistic") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_logist_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "log-logistic") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_llogist_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "probit") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_probit_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "log-probit") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_lprobit_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "multistage") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_multistage_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "qlinear") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_qlinear_f(fit_loop$parameters, test_doses)*weight
        }
        if (fit_loop$model == "weibull") {
          weight <- A$posterior_probs[ii]
          temp_f[ii, ] <- .dich_weibull_f(fit_loop$parameters, test_doses)*weight
        }
      }

      me <- colSums(temp_f)

      # Fitting line is not from sample- Need to double check with Matt

      lq <- apply(temp_f, 2, quantile, probs = qprob)
      uq <- apply(temp_f, 2, quantile, probs = 1 - qprob)

      temp_fit <- splinefun(test_doses, me)
      # Baseplot with minimal and maixmal dose with error bar
      plot_gg <- ggplot() +
        xlim(-max(doses) * 5, max(doses) * 5) +
        geom_errorbar(aes(x = doses, ymin = lerror, ymax = uerror), color = "grey") +
        ylim(c(min(Response, me, lq, uq) * 0.95, max(Response, me, lq, uq) * 1.05)) +
        labs(x = "Dose", y = "Proportion", title = "Model : Dichotomous MA, Fit type : Laplace") +
        theme_minimal()

      plot_gg <- plot_gg + geom_line(aes(x = test_doses, y = me), col = "blue", linewidth = 1.2) +
        geom_point(aes(x = doses, y = probs))


      # Laplace output doesn't have ..

      plot_gg <- plot_gg +
        geom_segment(aes(
          x = A$bmd[2], y = temp_fit(A$bmd[1]), xend = A$bmd[3],
          yend = temp_fit(A$bmd[1])
        ), color = "darkslategrey", linewidth = 1.2, alpha = 0.9)
      plot_gg <- plot_gg +
        annotate(
          geom = "text", x = A$bmd[2], y = temp_fit(A$bmd[1]),
          label = "[", size = 10, color = "darkslategrey", alpha = 0.9
        ) +
        annotate(
          geom = "text", x = A$bmd[3], y = temp_fit(A$bmd[1]),
          label = "]", size = 10, color = "darkslategrey", alpha = 0.9
        ) +
        annotate(
          geom = "point", x = A$bmd[1], y = temp_fit(A$bmd[1]),
          size = 5, color = "darkslategrey", shape = 17, alpha = 0.9
        )
      df <- NULL

      for (ii in 1:length(fit_idx)) {
        if (!is.finite(A$posterior_probs[ii])){
           A$posterior_probs[ii] = 0
        }
        if (A$posterior_probs[ii] > 0.05) {
          fit <- A[[ii]]
          if (fit$model == "hill") {
            f <- .dich_hill_f(fit$parameters, test_doses)
          }
          if (fit$model == "gamma") {
            f <- .dich_gamma_f(fit$parameters, test_doses)
          }
          if (fit$model == "logistic") {
            f <- .dich_logist_f(fit$parameters, test_doses)
          }
          if (fit$model == "log-logistic") {
            f <- .dich_llogist_f(fit$parameters, test_doses)
          }
          if (fit$model == "probit") {
            f <- .dich_probit_f(fit$parameters, test_doses)
          }
          if (fit$model == "log-probit") {
            f <- .dich_lprobit_f(fit$parameters, test_doses)
          }
          if (fit$model == "multistage") {
            f <- .dich_multistage_f(fit$parameters, test_doses)
          }
          if (fit$model == "qlinear") {
            f <- .dich_qlinear_f(fit$parameters, test_doses)
          }
          if (fit$model == "weibull") {
            f <- .dich_weibull_f(fit$parameters, test_doses)
          }

          col <- "coral3"
          temp_df <- data.frame(x_axis = test_doses, y_axis = f, cols = col, model_no = ii, alpha_lev = unname(A$posterior_probs[ii]))
          df <- rbind(df, temp_df)

          # SL Updated 06/18/21 -- Transparency update based on posterior probability and Y scale for dichotomous case
          temp_data <- df %>%
            filter(model_no == ii)

          plot_gg <- plot_gg +
            geom_line(data = temp_data, aes(x = x_axis, y = y_axis, color = cols), alpha = unique(temp_data$alpha_lev), show.legend = F) +
            theme_minimal()
        }
      }


      return(plot_gg + coord_cartesian(xlim = c(min(doses), max(doses)), expand = F))
    }
  }
}
