## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----load_data----------------------------------------------------------------
cont_data <- matrix(0, nrow = 5, ncol = 4)
colnames(cont_data) <- c("Dose", "Mean", "N", "SD")
cont_data[, 1] <- c(0, 50, 100, 200, 400)
cont_data[, 2] <- c(5.26, 5.76, 7.13, 9.24, 9.23)
cont_data[, 3] <- c(20, 20, 20, 20, 20)
cont_data[, 4] <- c(2.23, 1.47, 2.47, 2.24, 1.56)
Y <- cont_data[, 2:4]

## ----run_laplace_hill---------------------------------------------------------
library(ToxicR)
library(ggplot2)
hill_fit <- single_continuous_fit(cont_data[, "Dose"], Y,
  model_type = "hill"
)

## ----run_laplace_hill2--------------------------------------------------------
hill_fit$full_model
hill_fit$prior

## ----run_laplace_hill3, fig.height = 5, fig.width = 6-------------------------
hill_fit <- single_continuous_fit(cont_data[, "Dose"],
  cbind(cont_data[, "Mean"], cont_data[, "N"], cont_data[, "SD"]),
  model_type = "hill", distribution = "normal",
  fit_type = "mcmc"
)
hill_fit$full_model
plot(hill_fit)

## ----run_laplace_exp5, fig.height = 5, fig.width = 6--------------------------
exp5_fit <- single_continuous_fit(cont_data[, "Dose"], Y,
  model_type = "exp-5", distribution = "lognormal", fit_type = "laplace"
)
exp5_fit$full_model
plot(exp5_fit)

## ----run_laplace_hillad, fig.height = 5, fig.width = 6------------------------
hill_sd_fit <- single_continuous_fit(cont_data[, "Dose"], Y,
  model_type = "hill", distribution = "normal-ncv", fit_type = "mcmc",
  BMD_TYPE = "abs", BMR = 2
)
hill_sd_fit$full_model
plot(hill_sd_fit)

## ----run_laplace_hillsd, fig.height = 5, fig.width = 6------------------------
hill_sd_fit <- single_continuous_fit(cont_data[, "Dose"],
  cbind(cont_data[, "Mean"], cont_data[, "N"], cont_data[, "SD"]),
  model_type = "hill", distribution = "normal-ncv", fit_type = "laplace",
  BMD_TYPE = "sd", BMR = 1.5
)
hill_sd_fit$full_model
plot(hill_sd_fit)

## ----run_laplace_hillhybrid, fig.height = 5, fig.width = 6--------------------
hill_hybrid_fit <- single_continuous_fit(cont_data[, "Dose"],
  cbind(cont_data[, "Mean"], cont_data[, "N"], cont_data[, "SD"]),
  model_type = "hill", distribution = "normal-ncv", fit_type = "mcmc",
  BMD_TYPE = "hybrid", point_p = 0.025, BMR = 0.1
)
hill_hybrid_fit$full_model
plot(hill_hybrid_fit)

## ----run_laplace_hillrd, fig.height = 5, fig.width = 6------------------------
hill_rd_fit <- single_continuous_fit(cont_data[, "Dose"],
  cbind(cont_data[, "Mean"], cont_data[, "N"], cont_data[, "SD"]),
  model_type = "hill", distribution = "normal-ncv", fit_type = "mcmc",
  BMD_TYPE = "rel", BMR = 0.1, samples = 50000
)
hill_rd_fit$full_model
plot(hill_rd_fit)

## ----run_laplace_cprior-------------------------------------------------------
hill_sd_fit$prior

prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(0, 1, 0, 100),
  lnormprior(log(1), 0.4215, 0, 18),
  lnormprior(0, 1, 0, 100),
  normprior(0, 10, -100, 100)
)
p_hill_ncv <- create_continuous_prior(prior, "hill", "normal-ncv")

prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(0, 1, 0, 100),
  lnormprior(log(1), 0.4215, 0, 18),
  normprior(0, 10, -100, 100)
)
p_hill_norm <- create_continuous_prior(prior, "hill", "normal")

## ----run_laplace_cprior2, fig.height = 5, fig.width = 6-----------------------

hill_sd_a_fit <- single_continuous_fit(cont_data[, "Dose"],
  cbind(cont_data[, "Mean"], cont_data[, "N"], cont_data[, "SD"]),
  prior = p_hill_ncv,
  fit_type = "laplace",
  BMD_TYPE = "sd", BMR = 1.5
)

hill_sd_b_fit <- single_continuous_fit(cont_data[, "Dose"],
  cbind(cont_data[, "Mean"], cont_data[, "N"], cont_data[, "SD"]),
  prior = p_hill_norm,
  fit_type = "laplace",
  BMD_TYPE = "sd", BMR = 1.5
)

library(ggpubr)
figure <- ggarrange(plot(hill_sd_a_fit) + ggtitle(""),
  plot(hill_sd_b_fit) + ggtitle(""),
  labels = c("Prior NCV", "Prior Normal"),
  ncol = 1, nrow = 2
)

figure

## ----run_laplace_polynomial, fig.height = 5, fig.width = 6--------------------


poly_sd <- single_continuous_fit(cont_data[, "Dose"], Y,
  distribution = "normal", model_type = "polynomial",
  degree = 4,
  fit_type = "laplace",
  BMD_TYPE = "sd", BMR = 0.5
)

plot(poly_sd)

## ----run_laplace_MA_2---------------------------------------------------------
prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(0, 1, 0, 100),
  lnormprior(log(1), 0.4215, 0, 18),
  lnormprior(0, 1, 0, 100),
  normprior(0, 10, -100, 100)
)
p_hill_ncv <- create_continuous_prior(prior, "hill", "normal-ncv")

prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(0, 1, 0, 100),
  lnormprior(log(1), 0.4215, 0, 18),
  normprior(0, 10, -100, 100)
)
p_hill_norm <- create_continuous_prior(prior, "hill", "normal")

prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(0, 1, 0, 100),
  lnormprior(log(1), 0.4215, 0, 18),
  normprior(0, 10, -100, 100)
)
p_exp5_norm <- create_continuous_prior(prior, "exp-5", "normal")

prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(log(1), 0.4215, 0, 18),
  normprior(0, 10, -100, 100)
)
p_power_norm <- create_continuous_prior(prior, "power", "normal")

prior <- create_prior_list(
  normprior(0, 1, -100, 100),
  normprior(0, 1, -1e4, 1e4),
  lnormprior(log(1), 0.4215, 0, 18),
  normprior(0, 10, -100, 100)
)
p_exp3_norm <- create_continuous_prior(prior, "exp-3", "normal")
prior_list <- list(p_exp3_norm, p_hill_norm, p_exp5_norm, p_power_norm)

## ----run_laplace_MA_3, fig.height = 5, fig.width = 6--------------------------
ma_sd_mcmc_2 <- ma_continuous_fit(cont_data[, "Dose"], Y,
  fit_type = "laplace",
  BMD_TYPE = "sd", BMR = 0.5, samples = 50000, model_list = prior_list
)
plot(ma_sd_mcmc_2)

