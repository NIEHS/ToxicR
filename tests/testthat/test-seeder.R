context("GSL Seed")

test_that("Setting new seed for single continuous model", {
  set.seed(5981)
  M2 <- build_single_continuous_dataset()
  c <- single_continuous_fit(M2[, 1, drop = F], M2[, 2:4],
          BMR_TYPE = "sd", BMR = 1, ewald = T,
          distribution = "normal", fit_type = "laplace", model_type = "polynomial", degree = 4
     )
  validate_model(c, "Model: Polynomial Distribution: Normal", c(6.228232392127818e+00, -6.765515287902994e-02, -8.641725828795919e-05, 3.318759851015697e-06, -9.396333533806993e-09, -1.756114196384799e-02), c(14.52627182006836, 11.46386576089348, 18.40675539924204), tolerance = 10e-1)

  c <- single_continuous_fit(M2[, 1, drop = F], M2[, 2:4],
          BMR_TYPE = "sd", BMR = 1, ewald = T,
          distribution = "normal", fit_type = "laplace", model_type = "polynomial", degree = 4, seed = 11111
     )
  validate_model(c, "Model: Polynomial Distribution: Normal", c(6.228232351333556e+00, -6.765515484751268e-02, -8.641703828230548e-05, 3.318756799972366e-06, -9.396323504848345e-09, -1.756116168616254e-02), c(14.52627182006836, 11.45480350088830, 18.42131756988659), tolerance = 10e-1)
})