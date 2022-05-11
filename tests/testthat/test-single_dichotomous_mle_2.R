context("Single Dichotomous Models MLE")

test_that("Hill Laplace", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset_2()
     mData <- build_single_dichotomous_dataset_2()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],
                                model_type = "hill", fit_type = "laplace")
     validate_model2(c, "Model:  Hill", c(-1.673544, 3.459392,  5.833098, 1.462155), 
                    c(0.004217370,0.001890691,0.007291542),c(6.54215815,0.01053297))
})

test_that("Weibull Laplace", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset_2()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],
                              model_type = "weibull", fit_type = "laplace")
     validate_model2(c, "Model: Weibull", c(-1.7686897,0.7951624,14.4332211), 
                    c(0.0020553904,0.0007833962,0.0049656820 ),c(3.187 ,0.076))
})

test_that("Weibull MLE", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset_2()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],
                              model_type = "weibull", fit_type = "mle")
     validate_model2(c, "Model: Weibull", c(-1.8014517,0.7668152,13.6038641), 
                    c(0.0017663444,0.0006435623,0.0044451057),c(2.986,0.084))
})


test_that("Gamma MLE", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset_2()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],
                              model_type = "gamma", fit_type = "mle")
     validate_model2(c, "Model:  Gamma", c(-1.671756,1.000000, 24.832283), 
                    c(0.0020553904,0.0007833962,0.0049656820 ),c( 5.489 , 0.019))
})

test_that("Gamma MLE", {
     set.seed(5981)
     mData <- build_single_dichotomous_dataset_2()
     c = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],
                              model_type = "gamma", fit_type = "laplace")
     validate_model2(c, "Model:  Gamma",c( -1.7531361,0.7512065,18.8926902), 
                    c(0.004242885,0.003526430,0.006277881 ),c(3.125,0.085))
})

