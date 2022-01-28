context("Single Dichotomous Models MLE")

test_that("Hill", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mle")
     validate_model(c, "Model:  Hill", c(-3.15, -0.72, -8.184, 7.489), c(2.67, 1.34, 3.44))
})

test_that("Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit", fit_type = "mle")
     validate_model(c, "Model:  Probit", c(-1.22, 0.028), c(13.32, 10.99, 17.459))
})

test_that("Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit", fit_type = "mle")
     validate_model(c, "Model:  Log-Probit", c(-3.359, -1.55, 0.34), c(2.20, 0.578, 5.66))
})

test_that("Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull", fit_type = "mle")
     validate_model(c, "Model: Weibull", c(-3.367, 0.519, 0.04), c(2.02, 0.39, 5.989))
})