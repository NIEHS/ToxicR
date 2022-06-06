context("MA Continuous Models")


##
# Change 6/5/2022: Changed E. Wimberly's test conditions. No longer having the same
#                  test condition for all MA code. 
##
test_that("Laplace", {
     set.seed(5981)
     data <- build_ma_dataset()
     y = data[["y"]]
     doses = data[["doses"]]
     model_list <- build_model_list(y)
     AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y), 
                             fit_type = "laplace", BMD_TYPE = 'sd', BMR = 1)
     expect_equal(13, length(AA))
     expect_equal(setNames(c(8.568833,7.581585,9.725617), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance=10e-2)
     expect_equal(c(  2.490478e-04,9.997242e-01,4.919926e-32,1.029411e-31,4.184972e-30,
     0.068121e-08,2.383774e-05,2.865501e-06, 8.135752e-30,6.035466e-29), AA$posterior_probs, tolerance=10e-2)
     #generate_validation_code(AA)
        validate_model( AA$Individual_Model_1 ,  "Model: Hill Distribution: Normal" ,  c(10.5364575359962, 10.1321968130315, 25.9553804649865, 2.9312534185233, -1.78360527954745) ,  c(BMD = 8.81292036977129, BMDL = 7.77360319250264, BMDU = 10.038830163013) )
        validate_model( AA$Individual_Model_2 ,  "Model: Hill Distribution: Normal-NCV" ,  c(10.5387733317716, 10.1246980179567, 25.949625484858, 2.94190356756934, 0.217860255395499, -2.47791318851729) ,  c(BMD = 8.56879698769221, BMDL = 7.58530796701641, BMDU = 9.72728379864331) )
        validate_model( AA$Individual_Model_3 ,  "Model: Exponential-3 Distribution: Normal" ,  c(9.88000042605503, 0.00606034009825435, 0.501143334846178, 0.667522743525469) ,  c(BMD = 2.90979295969009, BMDL = 1.48524034218829, BMDU = 5.6128456958154) )
        validate_model( AA$Individual_Model_4 ,  "Model: Exponential-3 Distribution: Normal-NCV" ,  c(9.92091734166694, 0.00629062629297153, 0.526657833949637, 0.693441245848478, -1.1991037592188) ,  c(BMD = 2.64606848359108, BMDL = 1.38603949735297, BMDU = 5.00354043111249) )
        validate_model( AA$Individual_Model_5 ,  "Model: Exponential-3 Distribution: Log-Normal" ,  c(10.0536071161962, 0.00671019726121126, 0.594642014334531, -4.71379467965956) ,  c(BMD = 2.83083096146584, BMDL = 1.55282174831324, BMDU = 5.14804808012108) )
        validate_model( AA$Individual_Model_6 ,  "Model: Exponential-5 Distribution: Normal" ,  c(10.3572335733392, 0.0316363330323462, 0.678100646692081, 1.95623560832116, -1.52319151194146) ,  c(BMD = 6.66423290967941, BMDL = 5.69696817646146, BMDU = 7.88365455831351) )
        validate_model( AA$Individual_Model_7 ,  "Model: Exponential-5 Distribution: Normal-NCV" ,  c(10.4105947833681, 0.0318987775736513, 0.669930155029355, 2.03111716938036, 0.470804360632926, -2.88253758052939) ,  c(BMD = 6.60052299499512, BMDL = 5.48245053502279, BMDU = 7.62693632073755) )
        validate_model( AA$Individual_Model_8 ,  "Model: Exponential-5 Distribution: Log-Normal" ,  c(10.4453346497677, 0.0326706185835769, 0.661690178384068, 2.18753391383612, -6.9897057626201) ,  c(BMD = 6.47179484367371, BMDL = 5.45825730384399, BMDU = 7.81548945568862) )
        validate_model( AA$Individual_Model_9 ,  "Model: Power Distribution: Normal" ,  c(9.81455751339419, 0.607478195907743, 0.642833313623578, 0.581187397105843) ,  c(BMD = 3.41244720522803, BMDL = 1.98026031114979, BMDU = 5.88759957647574) )
        validate_model( AA$Individual_Model_10 ,  "Model: Power Distribution: Normal-NCV" ,  c(9.87814571011071, 0.528730524937364, 0.67558511920475, 0.687100269256969, -1.33060615543121) ,  c(BMD = 3.07466604865305, BMDL = 1.84384167291501, BMDU = 5.14769077561557) )

    
})


test_that("Vector Input", {
       
        data <- build_ma_dataset_2()
        y = data[["y"]]
        doses = data[["doses"]]
        model_list <- build_model_list(y)
        AA <- ma_continuous_fit(doses, y, model_list=model_list,
                                fit_type = "laplace", BMD_TYPE = 'sd', BMR = 1)
        expect_equal(13, length(AA))
        expect_equal(setNames(c(47.25482,37.39118,53.58764), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance=10e-2)
        expect_equal(c(  6.231843e-05, 2.980679e-02, 1.869943e-06, 2.638495e-06, 6.649184e-04,
                         1.160152e-04, 9.306290e-02, 8.762819e-01, 3.503513e-07, 2.615127e-07
        ), AA$posterior_probs, tolerance=10e-2)
        validate_model( AA$Individual_Model_1 ,  "Model: Hill Distribution: Normal" ,  c( 483.799766,-254.390813 ,70.933924,3.278107 ,7.587217) ,  c(BMD =44.16288, BMDL= 35.50511,BMDU= 52.29493) )
        validate_model( AA$Individual_Model_2 ,  "Model: Hill Distribution: Normal-NCV" ,  c( 476.193010, -234.611421  , 40.333183, 1.488576  ,1.733539 ,-3.101292) ,  c(BMD = 44.16288, BMDL = 35.50511, BMDU = 52.29493) )
        validate_model( AA$Individual_Model_3 ,  "Model: Exponential-3 Distribution: Normal" ,  c(4.843622e+02, 3.505894e-03, 6.737080e-01, 7.341241e+00) ,  c(BMD = 31.52798,BMDL=22.36762,BMDU=42.34427) )
        validate_model( AA$Individual_Model_4 ,  "Model: Exponential-3 Distribution: Normal-NCV" ,  c(  497.21082698, 0.00503978 , 1.13893843,1.78661758 ,-3.01411753) ,  c(BMD = 31.13267,22.17759,42.14889) )
        validate_model( AA$Individual_Model_5 ,  "Model: Exponential-3 Distribution: Log-Normal" ,  c(495.744954457 ,0.004958595,1.097434564 ,-4.305824373) ,  c(BMD = 31.73811,19.90481,38.76919) )
        validate_model( AA$Individual_Model_6 ,  "Model: Exponential-5 Distribution: Normal" ,  c( 483.72653627,0.0126337,-0.65566043,2.73463884,7.57828210) ,  c(BMD = 44.79610,35.74702,53.56469) )
        validate_model( AA$Individual_Model_7 ,  "Model: Exponential-5 Distribution: Normal-NCV" ,  c( 484.20344507 ,0.01286994,-0.64752041,2.92482439,1.80447278, -3.43029029) ,  c(BMD = 47.08060, BMDL = 39.29308, BMDU = 54.79945) )
        validate_model( AA$Individual_Model_8 ,  "Model: Exponential-5 Distribution: Log-Normal" ,  c( 481.25534823 ,0.01285451,-0.64296887,2.77758298,-4.60431736) ,  c(BMD = 47.28372, BMDL=37.22611, BMDU=53.38844) )
        validate_model( AA$Individual_Model_9 ,  "Model: Power Distribution: Normal" ,  c(494.907760,-1.559144,1.024907,7.826761) ,  c(BMD =  29.51604, BMDL = 19.61772, BMDU = 41.89054) )
        validate_model( AA$Individual_Model_10 ,  "Model: Power Distribution: Normal-NCV" ,  c( 498.8842946 , -2.7424269  , 0.9074885 ,1.7904453 , -2.9760857) ,  c(BMD = 29.28226, BMDL =19.38554,BMDU=42.00459) )
        
})


test_that("Plots", {
        set.seed(5981)
        data <- build_ma_dataset()
        y = data[["y"]]
        doses = data[["doses"]]
        model_list <- build_model_list(y)
        
        AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y), model_list=model_list,
                                fit_type = "laplace", BMD_TYPE = 'sd', BMR = 1)
        laplace_plot <- plot(AA)
        expect_identical(laplace_plot$labels$x, "Dose")
        expect_identical(laplace_plot$labels$y, "Response")
        expect_identical(laplace_plot$labels$title, "Continous MA fitting")
        
        laplace_cleveland <- cleveland_plot(AA)
        expect_identical(laplace_cleveland$labels$x, "Dose Level")
        expect_identical(laplace_cleveland$labels$title, "BMD Estimates by Each Model (Sorted by Posterior Probability)")
        
        AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                                fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1)
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
