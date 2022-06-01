context("MA Continuous Models")

test_that("No Model List", {
        set.seed(5981)
        data <- build_ma_dataset()
        y = data[["y"]]
        doses = data[["doses"]]
        AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y),
                                fit_type = "laplace", BMD_TYPE = 'sd', BMR = 1)
        expect_equal(13, length(AA))
        expect_equal(setNames(c(26.1, 13.1, 46.4), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance=10e-2)
        expect_equal(c(0.020170794, 0.088721617, 0.007274817, 0.010060694, 0.429308134, 0.032407617, 0.068653807, 0.330598226, 0.010474859, 0.002329435), AA$posterior_probs, tolerance=10e-2)
        #generate_validation_code(AA)
        validate_model( AA$Individual_Model_1 ,  "Model: Hill Distribution: Normal" ,  c(459.967845943872, -161.008754551028, 27.5188401639968, 1.65424081942613, 7.98247151057799) ,  c(BMD = 18.2375333442782, BMDL = 10.8444280330944, BMDU = 31.4585711838948) )
        validate_model( AA$Individual_Model_2 ,  "Model: Hill Distribution: Normal-NCV" ,  c(460.91986362233, -170.900638901152, 29.993036710078, 1.62508594445158, 1.77773636675909, -2.7526204573563) ,  c(BMD = 20.1864076616609, BMDL = 12.2303671904536, BMDU = 32.7000214684699) )
        validate_model( AA$Individual_Model_3 ,  "Model: Exponential-3 Distribution: Normal" ,  c(452.251132237139, 0.00360797504540128, 0.88089759918846, 8.01713036177101) ,  c(BMD = 27.3068406619132, BMDL = 13.8392818225457, BMDU = 48.1885256719032) )
        validate_model( AA$Individual_Model_4 ,  "Model: Exponential-3 Distribution: Normal-NCV" ,  c(453.72230059148, 0.00354755225529725, 0.858437854974876, 1.75670593800607, -2.57252569846925) ,  c(BMD = 28.7097525782883, BMDL = 15.0759417824207, BMDU = 48.2326956853754) )
        validate_model( AA$Individual_Model_5 ,  "Model: Exponential-3 Distribution: Log-Normal" ,  c(447.444965882228, 0.00323898710616777, 0.834695165188957, -3.99674729819358) ,  c(BMD = 33.555762283504, BMDL = 14.2225167196148, BMDU = 48.7244905437359) )
        validate_model( AA$Individual_Model_6 ,  "Model: Exponential-5 Distribution: Normal" ,  c(452.846141184847, 0.0118088430531916, -0.647844663710844, 1.03792464843229, 8.016117534184) ,  c(BMD = 26.0554229840636, BMDL = 14.0242097619719, BMDU = 45.6123344521286) )
        validate_model( AA$Individual_Model_7 ,  "Model: Exponential-5 Distribution: Normal-NCV" ,  c(453.642251519102, 0.0139479436503198, -0.57553385098482, 1.08215857051091, 1.75451155533531, -2.58703430717675) ,  c(BMD = 27.2294450551271, BMDL = 15.1763420518202, BMDU = 27.7930945677683) )
        validate_model( AA$Individual_Model_8 ,  "Model: Exponential-5 Distribution: Log-Normal" ,  c(447.941268941248, 0.0117667267210228, -0.624749228339207, 0.996583083830026, -4.01325233553182) ,  c(BMD = 31.3936345279217, BMDL = 14.3255087347062, BMDU = 45.7276526943185) )
        validate_model( AA$Individual_Model_9 ,  "Model: Power Distribution: Normal" ,  c(449.719459940181, -3.01167900801319, 0.850410029789607, 8.04702514307978) ,  c(BMD = 31.0261725414404, BMDL = 15.718578354615, BMDU = 53.5077103803206) )
        validate_model( AA$Individual_Model_10 ,  "Model: Power Distribution: Normal-NCV" ,  c(450.8323439329, -3.43254002670521, 0.822051300900347, 1.75784546770694, -2.56320443025547) ,  c(BMD = 32.2840716585567, BMDL = 17.0206144965613, BMDU = 53.1114395962907) )
        
})


test_that("Laplace", {
     set.seed(5981)
     data <- build_ma_dataset()
     y = data[["y"]]
     doses = data[["doses"]]
     model_list <- build_model_list(y)
     AA <- ma_continuous_fit(as.matrix(doses), as.matrix(y), model_list=model_list,
                             fit_type = "laplace", BMD_TYPE = 'sd', BMR = 1)
     expect_equal(13, length(AA))
     expect_equal(setNames(c(26.1, 13.1, 46.4), c("BMD", "BMDL", "BMDU")), AA$bmd, tolerance=10e-2)
     expect_equal(c(0.020170794, 0.088721617, 0.007274817, 0.010060694, 0.429308134, 0.032407617, 0.068653807, 0.330598226, 0.010474859, 0.002329435), AA$posterior_probs, tolerance=10e-2)
     #generate_validation_code(AA)
     validate_model( AA$Individual_Model_1 ,  "Model: Hill Distribution: Normal" ,  c(459.967845943872, -161.008754551028, 27.5188401639968, 1.65424081942613, 7.98247151057799) ,  c(BMD = 18.2375333442782, BMDL = 10.8444280330944, BMDU = 31.4585711838948) )
     validate_model( AA$Individual_Model_2 ,  "Model: Hill Distribution: Normal-NCV" ,  c(460.91986362233, -170.900638901152, 29.993036710078, 1.62508594445158, 1.77773636675909, -2.7526204573563) ,  c(BMD = 20.1864076616609, BMDL = 12.2303671904536, BMDU = 32.7000214684699) )
     validate_model( AA$Individual_Model_3 ,  "Model: Exponential-3 Distribution: Normal" ,  c(452.251132237139, 0.00360797504540128, 0.88089759918846, 8.01713036177101) ,  c(BMD = 27.3068406619132, BMDL = 13.8392818225457, BMDU = 48.1885256719032) )
     validate_model( AA$Individual_Model_4 ,  "Model: Exponential-3 Distribution: Normal-NCV" ,  c(453.72230059148, 0.00354755225529725, 0.858437854974876, 1.75670593800607, -2.57252569846925) ,  c(BMD = 28.7097525782883, BMDL = 15.0759417824207, BMDU = 48.2326956853754) )
     validate_model( AA$Individual_Model_5 ,  "Model: Exponential-3 Distribution: Log-Normal" ,  c(447.444965882228, 0.00323898710616777, 0.834695165188957, -3.99674729819358) ,  c(BMD = 33.555762283504, BMDL = 14.2225167196148, BMDU = 48.7244905437359) )
     validate_model( AA$Individual_Model_6 ,  "Model: Exponential-5 Distribution: Normal" ,  c(452.846141184847, 0.0118088430531916, -0.647844663710844, 1.03792464843229, 8.016117534184) ,  c(BMD = 26.0554229840636, BMDL = 14.0242097619719, BMDU = 45.6123344521286) )
     validate_model( AA$Individual_Model_7 ,  "Model: Exponential-5 Distribution: Normal-NCV" ,  c(453.642251519102, 0.0139479436503198, -0.57553385098482, 1.08215857051091, 1.75451155533531, -2.58703430717675) ,  c(BMD = 27.2294450551271, BMDL = 15.1763420518202, BMDU = 27.7930945677683) )
     validate_model( AA$Individual_Model_8 ,  "Model: Exponential-5 Distribution: Log-Normal" ,  c(447.941268941248, 0.0117667267210228, -0.624749228339207, 0.996583083830026, -4.01325233553182) ,  c(BMD = 31.3936345279217, BMDL = 14.3255087347062, BMDU = 45.7276526943185) )
     validate_model( AA$Individual_Model_9 ,  "Model: Power Distribution: Normal" ,  c(449.719459940181, -3.01167900801319, 0.850410029789607, 8.04702514307978) ,  c(BMD = 31.0261725414404, BMDL = 15.718578354615, BMDU = 53.5077103803206) )
     validate_model( AA$Individual_Model_10 ,  "Model: Power Distribution: Normal-NCV" ,  c(450.8323439329, -3.43254002670521, 0.822051300900347, 1.75784546770694, -2.56320443025547) ,  c(BMD = 32.2840716585567, BMDL = 17.0206144965613, BMDU = 53.1114395962907) )
     
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
