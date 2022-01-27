context("Single Continuous Models Laplace")

test_that("Normal Ewald Hill", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="laplace",model_type = "hill",degree = 4)
     validate_model(c, "Model: Hill Distribution: Normal", c(6.15, -5.33, 39.89, 3.12, -0.15), c(24.01, 20.12, 28.66))
})

test_that("Normal Ewald exp-3", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="laplace",model_type = "exp-3",degree = 4)
     validate_model(c, "Model: Exponential-3 Distribution: Normal", c(6.149, 0.0159, 1.414, 0.049), c(18.86, 13.949, 25.51))
})

test_that("Normal Ewald exp-5", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="laplace",model_type = "exp-5",degree = 4)
     validate_model(c, "Model: Exponential-5 Distribution: Normal", c(5.99, 0.02, -1.89, 2.487, -0.17), c(25.2, 20.3, 31.4))
})

test_that("Normal Ewald power", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="laplace",model_type = "power",degree = 4)
     validate_model(c, "Model: Power Distribution: Normal", c(6.077, -0.380, 0.517, 0.339), c(8.98, 4.71, 17.116))
})

test_that("Normal Ewald polynomial", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="laplace",model_type = "polynomial",degree = 4)
     validate_model(c, "Model: Polynomial Distribution: Normal", c(6.23, -0.062, -0.000085, 0.0000033, -0.00000000936, -0.013), c(14.6, 11.5, 18.5))
})