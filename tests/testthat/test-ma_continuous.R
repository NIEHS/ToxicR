context("MA Continuous Models")


##
# Change 6/5/2022: Changed E. Wimberly's test conditions. No longer having the same
#                  test condition for all MA code.
##
test_that("Laplace", {
        set.seed(5981)
        data <- build_ma_dataset()
        y <- data[["y"]]
        doses <- data[["doses"]]
        model_list <- build_model_list(y)
        AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y),
                fit_type = "laplace", BMR_TYPE = "sd", BMR = 1,EFSA=FALSE
        )
        expect_equal(13, length(AA))
        expect_equal(setNames(c(6.432882, 5.456776, 7.819051), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance = 10e-2)
        expect_equal(setNames(c(
           1.371068e-12, 2.343204e-04, 9.652573e-29, 5.895678e-32, 1.480816e-24, 5.629987e-16, 1.937811e-09,
           9.997657e-01, 2.563253e-33, 1.819908e-29),
           c('hill_normal', 'hill_normal-ncv', 'exp-3_normal', 'exp-3_normal-ncv', 'exp-3_lognormal', 'exp-5_normal', 'exp-5_normal-ncv', 'exp-5_lognormal', 'power_normal', 'power_normal-ncv')),
           AA$posterior_probs, tolerance = 10e-3)
        # generate_validation_code(AA)
        validate_model( AA$Indiv_hill_normal ,  "Model: Hill Distribution: Normal" ,  c(10.5336679441343, 10.1406923686382, 25.9628913569909, 2.92082302495687, -1.58762629617883) ,  c(BMD = 9.09194176320145, BMDL = 7.94835310278431, BMDU = 10.453730877853) )
        validate_model( AA$`Indiv_hill_normal-ncv` ,  "Model: Hill Distribution: Normal-NCV" ,  c(10.5358381990029, 10.1336040160307, 25.9571953249711, 2.93024372184542, 0.096233868317774, -1.97358149131162) ,  c(BMD = 8.8674227836588, BMDL = 7.79233267848719, BMDU = 10.1379381094965) )
        validate_model( AA$`Indiv_exp-3_normal` ,  "Model: Exponential-3 Distribution: Normal" ,  c(9.88651394682159, 0.00605728233962403, 0.50215419205304, 0.708257513692966) ,  c(BMD = 3.04486677050591, BMDL = 1.54246823432622, BMDU = 5.91744535584554) )
        validate_model( AA$`Indiv_exp-3_normal-ncv` ,  "Model: Exponential-3 Distribution: Normal-NCV" ,  c(9.91444327139625, 0.00616455639933146, 0.516273721648633, 0.331008185094923, -0.0814104874597511) ,  c(BMD = 3.23001891374588, BMDL = 1.62523168417107, BMDU = 6.3323948292675) )
        validate_model( AA$`Indiv_exp-3_lognormal` ,  "Model: Exponential-3 Distribution: Log-Normal" ,  c(10.0536073351969, 0.00671019676358537, 0.594642026870305, -4.71379466264071) ,  c(BMD = 2.8308317065239, BMDL = 1.55279522712495, BMDU = 5.14804943505705) )
        validate_model( AA$`Indiv_exp-5_normal` ,  "Model: Exponential-5 Distribution: Normal" ,  c(10.3506333415941, 0.0315674536891819, 0.679299421765513, 1.94736736795438, -1.32829634007905) ,  c(BMD = 6.97416514158249, BMDL = 5.91422854869825, BMDU = 8.37825850036936) )
        validate_model( AA$`Indiv_exp-5_normal-ncv` ,  "Model: Exponential-5 Distribution: Normal-NCV" ,  c(10.375276768482, 0.0317126958629768, 0.675465222273146, 1.98258411488366, 0.272276041055903, -2.17702052717982) ,  c(BMD = 6.77031278610229, BMDL = 5.72202795631237, BMDU = 8.02395634143899) )
        validate_model( AA$`Indiv_exp-5_lognormal` ,  "Model: Exponential-5 Distribution: Log-Normal" ,  c(10.440056037444, 0.0326170917489223, 0.662682229270099, 2.17773617305393, -6.98862503700197) ,  c(BMD = 6.43275082111359, BMDL = 5.45804426103146, BMDU = 7.81594398797255) )
        validate_model( AA$Indiv_power_normal ,  "Model: Power Distribution: Normal" ,  c(9.83030261567371, 0.598353607312653, 0.645939852257835, 0.759828909905049) ,  c(BMD = 3.98778882856081, BMDL = 2.24179312948605, BMDU = 7.11222865500403) )
        validate_model( AA$`Indiv_power_normal-ncv` ,  "Model: Power Distribution: Normal-NCV" ,  c(9.84909259080531, 0.564614132759564, 0.659824365338973, 0.329109555698286, -0.220809149647623) ,  c(BMD = 3.55885856887591, BMDL = 2.05664907507987, BMDU = 6.18269998269575) )
})


test_that("Vector Input", {
        data <- build_ma_dataset_2()
        y <- data[["y"]]
        doses <- data[["doses"]]
        model_list <- build_model_list(y)
        AA <- ma_continuous_fit(doses, y,
                model_list = model_list,
                fit_type = "laplace", BMR_TYPE = "sd", BMR = 1,EFSA=FALSE
        )
        expect_equal(13, length(AA))
        expect_equal(setNames(c(47.25482, 37.39118, 53.58764), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance = 10e-2)
        expect_equal(setNames(c(
                6.231843e-05, 2.980679e-02, 1.869943e-06, 2.638495e-06, 6.649184e-04,
                1.160152e-04, 9.306290e-02, 8.762819e-01, 3.503513e-07, 2.615127e-07
        ),c('hill_normal', 'hill_normal-ncv', 'exp-3_normal', 'exp-3_normal-ncv', 'exp-3_lognormal', 'exp-5_normal', 'exp-5_normal-ncv', 'exp-5_lognormal', 'power_normal', 'power_normal-ncv')),
        AA$posterior_probs, tolerance = 10e-2)
        validate_model(AA$Indiv_hill_normal, "Model: Hill Distribution: Normal", c(483.799766, -254.390813, 70.933924, 3.278107, 7.587217), c(BMD = 44.16288, BMDL = 35.50511, BMDU = 52.29493))
        validate_model(AA$`Indiv_hill_normal-ncv`, "Model: Hill Distribution: Normal-NCV", c(476.193010, -234.611421, 40.333183, 1.488576, 1.733539, -3.101292), c(BMD = 44.16288, BMDL = 35.50511, BMDU = 52.29493))
        validate_model(AA$`Indiv_exp-3_normal`, "Model: Exponential-3 Distribution: Normal", c(4.843622e+02, 3.505894e-03, 6.737080e-01, 7.341241e+00), c(BMD = 31.52798, BMDL = 22.36762, BMDU = 42.34427))
        validate_model(AA$`Indiv_exp-3_normal-ncv`, "Model: Exponential-3 Distribution: Normal-NCV", c(497.21082698, 0.00503978, 1.13893843, 1.78661758, -3.01411753), c(BMD = 31.13267, 22.17759, 42.14889))
        validate_model(AA$`Indiv_exp-3_lognormal`, "Model: Exponential-3 Distribution: Log-Normal", c(495.744954457, 0.004958595, 1.097434564, -4.305824373), c(BMD = 31.73811, 19.90481, 38.76919))
        validate_model(AA$`Indiv_exp-5_normal`, "Model: Exponential-5 Distribution: Normal", c(483.72653627, 0.0126337, -0.65566043, 2.73463884, 7.57828210), c(BMD = 44.79610, 35.74702, 53.56469))
        validate_model(AA$`Indiv_exp-5_normal-ncv`, "Model: Exponential-5 Distribution: Normal-NCV", c(484.20344507, 0.01286994, -0.64752041, 2.92482439, 1.80447278, -3.43029029), c(BMD = 47.08060, BMDL = 39.29308, BMDU = 54.79945))
        validate_model(AA$`Indiv_exp-5_lognormal`, "Model: Exponential-5 Distribution: Log-Normal", c(481.25534823, 0.01285451, -0.64296887, 2.77758298, -4.60431736), c(BMD = 47.28372, BMDL = 37.22611, BMDU = 53.38844))
        validate_model(AA$Indiv_power_normal, "Model: Power Distribution: Normal", c(494.907760, -1.559144, 1.024907, 7.826761), c(BMD = 29.51604, BMDL = 19.61772, BMDU = 41.89054))
        validate_model(AA$`Indiv_power_normal-ncv`, "Model: Power Distribution: Normal-NCV", c(498.8842946, -2.7424269, 0.9074885, 1.7904453, -2.9760857), c(BMD = 29.28226, BMDL = 19.38554, BMDU = 42.00459))
})


test_that("Plots", {
        set.seed(5981)
        data <- build_ma_dataset()
        y <- data[["y"]]
        doses <- data[["doses"]]
        model_list <- build_model_list(y)

        AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y),
                model_list = model_list,
                fit_type = "laplace", BMR_TYPE = "sd", BMR = 1,EFSA=FALSE
        )
        laplace_plot <- plot(AA)
        expect_identical(laplace_plot$labels$x, "Dose")
        expect_identical(laplace_plot$labels$y, "Response")
        expect_identical(laplace_plot$labels$title, "Continous MA fitting")

        laplace_cleveland <- cleveland_plot(AA)
        expect_identical(laplace_cleveland$labels$x, "Dose Level")
        expect_identical(laplace_cleveland$labels$title, "BMD Estimates by Each Model (Sorted by Posterior Probability)")

        AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y),
                model_list = model_list,
                fit_type = "mcmc", BMR_TYPE = "sd", BMR = 1
        )
        mcmc_plot <- plot(AA)
        expect_identical(mcmc_plot$labels$x, "Dose")
        expect_identical(mcmc_plot$labels$y, "Response")
        expect_identical(mcmc_plot$labels$title, "Continous MA fitting")

        mcmc_cleveland <- cleveland_plot(AA)
        expect_identical(mcmc_cleveland$labels$x, "Dose Level")
        expect_identical(mcmc_cleveland$labels$title, "BMD Estimates by Each Model (Sorted by Posterior Probability)")

        mcmc_density <- MAdensity_plot(AA)
        expect_identical(mcmc_density$labels$x, "Dose Level (Dotted Line : MA BMD)")
        expect_identical(mcmc_density$labels$title, "Density plots for each fitted model (Fit type: MCMC)")
})
