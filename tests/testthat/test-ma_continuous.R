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
                fit_type = "laplace", BMR_TYPE = "sd", BMR = 1, EFSA = FALSE
        )
  expect_equal(13, length(AA))
  expect_equal(setNames(c(8.492525612073692, 5.836868675655248, 9.991342142246316), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance = 10e-2)
  expect_equal(setNames(c(
      4.105796851001509e-09, 7.016067040180166e-01, 2.890703321219598e-29, 1.748444162039573e-32, 4.433970447733191e-25, 1.607315649069202e-16, 4.874538760946350e-10,
      2.983932913887326e-01, 3.919246629344100e-33, 2.605436006966734e-29),
      c('hill_normal', 'hill_normal-ncv', 'exp-3_normal', 'exp-3_normal-ncv', 'exp-3_lognormal', 'exp-5_normal', 'exp-5_normal-ncv', 'exp-5_lognormal', 'power_normal', 'power_normal-ncv')),
      AA$posterior_probs, tolerance = 10e-3)
  # generate_validation_code(AA)
  validate_model(AA$`Indiv_hill_normal-ncv`, "Model: Hill Distribution: Normal-NCV", c(10.535838188644, 10.1336038451559, 25.9571944239826, 2.93024384912606, 0.0962323824892475, -1.97357757945706), c(BMD = 8.867423594475, BMDL = 7.792365015406, BMDU = 10.137939036486))
  validate_model(AA$Indiv_hill_normal, "Model: Hill Distribution: Normal", c(10.5336680740673, 10.1406918876008, 25.9628909121109, 2.92082355405902, -1.58762628346099), c(BMD = 9.091943510830, BMDL = 7.948330052503, BMDU = 10.453733757487))
  validate_model(AA$`Indiv_exp-3_normal`, "Model: Exponential-3 Distribution: Normal", c(9.88651312621504, 0.00605728353693025, 0.502154128307133, 0.708257514343457), c(BMD = 3.044865280390, BMDL = 1.542469608265, BMDU = 5.917305522635))
  validate_model(AA$`Indiv_exp-3_normal-ncv`, "Model: Exponential-3 Distribution: Normal-NCV", c(9.9144433272954, 0.00616455621043127, 0.516273717108163, 0.331008110199253, -0.0814103219567423), c(BMD = 3.230018913746, BMDL = 1.625237043886, BMDU = 6.332500836495))
  validate_model(AA$`Indiv_exp-3_lognormal`, "Model: Exponential-3 Distribution: Log-Normal", c(10.0536071579414, 0.00671019709770269, 0.594642014544935, -4.71379464978166), c(BMD = 2.830830961466, BMDL = 1.552794818438, BMDU = 5.148054709149))
  validate_model(AA$`Indiv_exp-5_normal`, "Model: Exponential-5 Distribution: Normal", c(10.3198495594875, 0.0315697493220312, 0.682793711047602, 1.92947484997103, -1.32354469345963), c(BMD = 6.870178878307, BMDL = 5.897494207110, BMDU = 7.012391581088))
  validate_model(AA$`Indiv_exp-5_normal-ncv`, "Model: Exponential-5 Distribution: Normal-NCV", c(10.3054183950297, 0.0316962480160407, 0.683696153196204, 1.93794411345487, 0.268202136937105, -2.15510619282155), c(BMD = 6.80, BMDL = 5.677321330143, BMDU = 8.10), tolerance = 3e-1)
  validate_model(AA$`Indiv_exp-5_lognormal`, "Model: Exponential-5 Distribution: Log-Normal", c(10.4408065945634, 0.0326044739570929, 0.66264718304348, 2.17637123560348, -6.98830954478905), c(BMD = 6.429639458656, BMDL = 5.457942006254, BMDU = 7.816159010457))
  validate_model(AA$Indiv_power_normal, "Model: Power Distribution: Normal", c(9.83030248873981, 0.598353667869175, 0.645939833044569, 0.759828898643816), c(BMD = 3.987788333070, BMDL = 2.241792850939, BMDU = 7.112227771295))
  validate_model(AA$`Indiv_power_normal-ncv`, "Model: Power Distribution: Normal-NCV", c(9.84909255672426, 0.56461411919398, 0.659824373649294, 0.329109984548888, -0.220810357935071), c(BMD = 3.558858025389, BMDL = 2.056648761001, BMDU = 6.182699038512))
})


test_that("Vector Input", {
  data <- build_ma_dataset_2()
  y <- data[["y"]]
  doses <- data[["doses"]]
  model_list <- build_model_list(y)
  AA <- ma_continuous_fit(doses, y,
                model_list = model_list,
                fit_type = "laplace", BMR_TYPE = "sd", BMR = 1, EFSA = FALSE
        )
  expect_equal(13, length(AA))
  expect_equal(setNames(c(47.25482, 37.39118, 53.58764), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance = 10e-2)
  expect_equal(setNames(c(
                6.231843e-05, 2.980679e-02, 1.869943e-06, 2.638495e-06, 6.649184e-04,
                1.160152e-04, 9.306290e-02, 8.762819e-01, 3.503513e-07, 2.615127e-07
        ), c('hill_normal', 'hill_normal-ncv', 'exp-3_normal', 'exp-3_normal-ncv', 'exp-3_lognormal', 'exp-5_normal', 'exp-5_normal-ncv', 'exp-5_lognormal', 'power_normal', 'power_normal-ncv')),
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
                fit_type = "laplace", BMR_TYPE = "sd", BMR = 1, EFSA = FALSE
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
