context("Single Dichotomous Models Laplace")

test_that("Hill", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "laplace")
     validate_model(c, "Model:  Hill", c(-3.15, -0.54, -2.215, 1.38), c(2.43, 1.19, 5.66))
})

test_that("Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit", fit_type = "laplace")
     validate_model(c, "Model:  Probit", c(-1.22, 0.028), c(13.32, 10.99, 17.459))
})

test_that("Log Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit", fit_type = "laplace")
     validate_model(c, "Model:  Log-Probit", c(-2.76, -2.03, 0.487), c(4.659, 1.88, 12.00))
})

test_that("Weibull", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull", fit_type = "laplace")
     validate_model(c, "Model: Weibull", c(-2.86, 0.66, 0.04), c(3.81, 1.26, 9.68))
})

test_that("Vector Input", {
        set.seed(5981)
        mData <- build_single_dichotomous_dataset()
        D <- as.double(mData[,1])
        dim(D) <- c(nrow(mData),1)
        Y <- as.double(mData[,2])
        dim(Y) <- c(nrow(mData),1)
        N <- as.double(mData[,3])
        dim(N) <- c(nrow(mData),1)
        c = single_dichotomous_fit(D, Y, N, model_type = "hill", fit_type = "laplace")
        validate_model(c, "Model:  Hill", c(-3.15, -0.54, -2.215, 1.38), c(2.43, 1.19, 5.66))
})

test_that("Plots", {
        set.seed(5981)
        mData <- build_single_dichotomous_dataset()
        c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "laplace")
        laplace_plot <- plot(c)
        expect_identical(laplace_plot$labels$x, "Dose")
        expect_identical(laplace_plot$labels$y, "Proportion")
        expect_identical(laplace_plot$labels$title, "Model:  Hill")
        
        c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull", fit_type = "laplace")
        laplace_plot <- plot(c)
        expect_identical(laplace_plot$labels$x, "Dose")
        expect_identical(laplace_plot$labels$y, "Proportion")
        expect_identical(laplace_plot$labels$title, "Model: Weibull")
})