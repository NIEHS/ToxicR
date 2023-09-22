context("MA Dichotomous Models")

test_that("Defaults", {
        set.seed(5981)
        mData <- build_single_dichotomous_dataset()

        AA <- ma_dichotomous_fit(mData[, 1], mData[, 2], mData[, 3])

        expect_equal(c("BMDdichotomous_MA", "BMDdichotomous_MA_laplace"), class(AA))
        expect_equal(13, length(AA))
        expect_equal(setNames(c(4.17, 1.33, 12.318), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance = 10e-2)
        expect_equal(setNames(c(0.427802462, 0.021635893, 0.022796957, 0.109884198, 0.006866267, 0.001184459, 0.025631532, 0.291857962, 0.092340269), 
                              c("hill_","gamma_","logistic_","log-logistic_","log-probit_","multistage_","probit_","qlinear_","weibull_")), 
                              AA$posterior_probs, tolerance = 10e-2)
        # generate_validation_code(AA)
        validate_model(AA$Indiv_hill_, "Model:  Hill", c(-3.1514129532799, -0.542504910549596, -2.21526647189887, 1.38193548558345), c(2.43689558842919, 1.19423317401394, 5.66302904194827))
        validate_model(AA$Indiv_gamma_, "Model:  Gamma", c(-2.58597525196282, 0.859030175600515, 0.0103251850555099), c(6.46707681632849, 2.71101742377789, 13.1872288945701))
        validate_model(AA$Indiv_logistic_, "Model:  Logistic", c(-2.0069607207327, 0.0463222123387509), c(14.2821369829316, 11.7611052273976, 18.9167440555661))
        validate_model(AA$`Indiv_log-logistic_`, "Model:  Log-Logistic", c(-2.90509719095083, -3.22382968699211, 0.770460283795502), c(3.79034439690537, 1.37987439309905, 9.29526839843571))
        validate_model(AA$`Indiv_log-probit_`, "Model:  Log-Probit", c(-2.76563341583274, -2.03205435911903, 0.487681754992638), c(4.65955128862954, 1.88253828334611, 12.0003951662748))
        validate_model(AA$Indiv_multistage_, "Model:  Multistage", c(-2.50756489338324, 0.0122508853947786, 8.85365951190481e-05), c(8.12334001064301, 6.09589901365889, 11.2063977029104))
        validate_model(AA$Indiv_probit_, "Model:  Probit", c(-1.22897249357605, 0.0286868893347652), c(13.3264085482876, 10.9912280971799, 17.4597665445259))
        validate_model(AA$Indiv_qlinear_, "Model:  Quantal-Linear", c(-2.45423133972787, 0.0132202483334582), c(7.9696321128232, 5.93259065501385, 11.5959106855121))
        validate_model(AA$Indiv_weibull_, "Model: Weibull", c(-2.86014800998383, 0.661116573139896, 0.0434668812243724), c(3.81611266939088, 1.26068620197916, 9.68226547559964))
})

test_that("Vector Inputs", {
        set.seed(5981)
        mData <- build_single_dichotomous_dataset()
        D <- as.double(mData[, 1])
        dim(D) <- c(nrow(mData), 1)
        Y <- as.double(mData[, 2])
        dim(Y) <- c(nrow(mData), 1)
        N <- as.double(mData[, 3])
        dim(N) <- c(nrow(mData), 1)
        AA <- ma_dichotomous_fit(D, Y, N)

        expect_equal(c("BMDdichotomous_MA", "BMDdichotomous_MA_laplace"), class(AA))
        expect_equal(13, length(AA))
        expect_equal(setNames(c(4.17, 1.33, 12.318), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance = 10e-2)
        expect_equal(setNames(c(0.427802462, 0.021635893, 0.022796957, 0.109884198, 0.006866267, 0.001184459, 0.025631532, 0.291857962, 0.092340269),
                              c("hill_","gamma_","logistic_","log-logistic_","log-probit_","multistage_","probit_","qlinear_","weibull_")), 
                     AA$posterior_probs, tolerance = 10e-2)
        # generate_validation_code(AA)
        validate_model(AA$Indiv_hill_, "Model:  Hill", c(-3.1514129532799, -0.542504910549596, -2.21526647189887, 1.38193548558345), c(2.43689558842919, 1.19423317401394, 5.66302904194827))
        validate_model(AA$Indiv_gamma_, "Model:  Gamma", c(-2.58597525196282, 0.859030175600515, 0.0103251850555099), c(6.46707681632849, 2.71101742377789, 13.1872288945701))
        validate_model(AA$Indiv_logistic_, "Model:  Logistic", c(-2.0069607207327, 0.0463222123387509), c(14.2821369829316, 11.7611052273976, 18.9167440555661))
        validate_model(AA$`Indiv_log-logistic_`, "Model:  Log-Logistic", c(-2.90509719095083, -3.22382968699211, 0.770460283795502), c(3.79034439690537, 1.37987439309905, 9.29526839843571))
        validate_model(AA$`Indiv_log-probit_`, "Model:  Log-Probit", c(-2.76563341583274, -2.03205435911903, 0.487681754992638), c(4.65955128862954, 1.88253828334611, 12.0003951662748))
        validate_model(AA$Indiv_multistage_, "Model:  Multistage", c(-2.50756489338324, 0.0122508853947786, 8.85365951190481e-05), c(8.12334001064301, 6.09589901365889, 11.2063977029104))
        validate_model(AA$Indiv_probit_, "Model:  Probit", c(-1.22897249357605, 0.0286868893347652), c(13.3264085482876, 10.9912280971799, 17.4597665445259))
        validate_model(AA$Indiv_qlinear_, "Model:  Quantal-Linear", c(-2.45423133972787, 0.0132202483334582), c(7.9696321128232, 5.93259065501385, 11.5959106855121))
        validate_model(AA$Indiv_weibull_, "Model: Weibull", c(-2.86014800998383, 0.661116573139896, 0.0434668812243724), c(3.81611266939088, 1.26068620197916, 9.68226547559964))
})

test_that("Plots", {
        set.seed(5981)
        mData <- build_single_dichotomous_dataset()
        AA <- ma_dichotomous_fit(mData[, 1], mData[, 2], mData[, 3])

        dichotomous_plot <- plot(AA)
        expect_identical(dichotomous_plot$labels$x, "Dose")
        expect_identical(dichotomous_plot$labels$y, "Proportion")
        expect_identical(dichotomous_plot$labels$title, "Model : Dichotomous MA, Fit type : Laplace")

        dichotomous_cleveland <- cleveland_plot(AA)
        expect_identical(dichotomous_cleveland$labels$x, "Dose Level")
        expect_identical(dichotomous_cleveland$labels$title, "BMD Estimates by Each Model (Sorted by Posterior Probability)")

        AA <- ma_dichotomous_fit(mData[, 1], mData[, 2], mData[, 3], fit_type = "mcmc")
        dichotomous_plot <- plot(AA)
        expect_identical(dichotomous_plot$labels$x, "Dose")
        expect_identical(dichotomous_plot$labels$y, "Proportion")
        # TODO should fit type MCMC be in the title?
        expect_identical(dichotomous_plot$labels$title, "Model : Dichotomous MA")

        dichotomous_cleveland <- cleveland_plot(AA)
        expect_identical(dichotomous_cleveland$labels$x, "Dose Level")
        expect_identical(dichotomous_cleveland$labels$title, "BMD Estimates by Each Model (Sorted by Posterior Probability)")

        density_plot <- MAdensity_plot(AA)
        expect_identical(density_plot$labels$x, "Dose Level (Dotted Line : MA BMD)")
        expect_identical(density_plot$labels$title, "Density plots for each fitted model (Fit type: MCMC)")
})
