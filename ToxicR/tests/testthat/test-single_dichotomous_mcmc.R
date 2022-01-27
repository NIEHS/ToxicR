context("Single Dichotomous Models MCMC")

test_that("Hill", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")
     validate_model(c,  "Model:  Hill" ,  c(-3.1514129532799, -0.542504910549596, 2.61668140290251, 1.38193548558345) ,  c(BMD = 3.22355133399257, BMDL = 1.32570638651708, BMDU = 7.2168846887097) )
})

test_that("Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit", fit_type = "mcmc")
     validate_model(c,  "Model:  Probit" ,  c(-1.22897249357605, 0.946667348047251) ,  c(BMD = 13.7570231637969, BMDL = 11.0421600688424, BMDU = 17.634673852233) )
})

test_that("Log Probit", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit", fit_type = "mcmc")
     validate_model(c,  "Model:  Log-Probit" ,  c(-2.76563341583274, -0.326871415198025, 0.487681754992638) ,  c(BMD = 7.03560338819198, BMDL = 2.26582400694548, BMDU = 15.4151199324048) )
})

test_that("Weibull", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull", fit_type = "mcmc")
     validate_model(c,  "Model: Weibull" ,  c(-2.86014800998383, 0.661116573139896, 0.438604630748536) ,  c(BMD = 5.1170096470376, BMDL = 1.46662721211377, BMDU = 11.8009898041251) )
})

test_that("Plots", {
        set.seed(5981)
        mData <- build_single_dichotomous_dataset()
        c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull", fit_type = "mcmc")
        mcmc_plot <- plot(c)
        expect_identical(mcmc_plot$labels$x, "Dose")
        expect_identical(mcmc_plot$labels$y, "Proportion")
        #TODO should the title have the distribution name?
        expect_identical(mcmc_plot$labels$title, "Model: Weibull,  Fit Type: MCMC")
        c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")
        mcmc_plot <- plot(c)
        expect_identical(mcmc_plot$labels$x, "Dose")
        expect_identical(mcmc_plot$labels$y, "Proportion")
        #TODO should the title have the distribution name?
        expect_identical(mcmc_plot$labels$title, "Model:  Hill,  Fit Type: MCMC")
})