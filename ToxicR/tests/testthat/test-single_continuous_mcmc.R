context("Single Continuous Models MCMC")

test_that("Normal Ewald Hill", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="mcmc",model_type = "hill",degree = 4)
     validate_model(c,  "Model: Hill Distribution: Normal" ,  c(6.0730964791569, -5.32930641539168, 39.8852409296944, 3.11573367222951, -0.188967044405645) ,  c(BMD = 24.527770998594, BMDL = 20.2417034249427, BMDU = 29.2708718586858) )
})

test_that("Normal Ewald exp-3", {
        #set.seed(5981)
        set.seed(5983)
        M2 <- build_single_continuous_dataset()
        c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "exp-3",degree = 4)
        validate_model(c,  "Model: Exponential-3 Distribution: Normal" ,  c(6.14912312280982, 0.0159199199226095, 1.4148159969476, 0.0493469100167165) ,  c(BMD = 19.2543569803238, BMDL = 13.9079211711884, BMDU = 25.7881819248199) )
})

test_that("Normal Ewald exp-5", {
        set.seed(5981)
        M2 <- build_single_continuous_dataset()
        c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "exp-5",degree = 4)
        validate_model(c,  "Model: Exponential-5 Distribution: Normal" ,  c(5.99508753796603, 0.0206933938641134, -1.89230053075426, 2.48798682208172, -0.171452987565699) ,  c(BMD = 25.3909630537033, BMDL = 19.9345767259598, BMDU = 31.0367735505104) )
})

test_that("Normal Ewald power", {
        set.seed(5981)
        M2 <- build_single_continuous_dataset()
        c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "power",degree = 4)
        validate_model(c,  "Model: Power Distribution: Normal" ,  c(6.07789287920542, -0.38081482983838, 0.517160540931021, 0.339528184806466) ,  c(BMD = 9.54618534364222, BMDL = 4.76945861911754, BMDU = 18.141129083891) )
})

test_that("Normal Ewald polynomial", {
        set.seed(5981)
        M2 <- build_single_continuous_dataset()
        c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "polynomial",degree = 4)
        validate_model(c,  "Model: Polynomial Distribution: Normal" ,  c(6.22823399110153, -0.0676790046609912, -8.53912855307723e-05, 3.30743840069754e-06, -9.36239323304629e-09, -0.0134111696312258) ,  c(BMD = 14.9544341564178, BMDL = 11.5303337335587, BMDU = 20.5506707191467) )
})

test_that("Vector Inputs", {
        set.seed(5981)
        M2 <- build_single_continuous_dataset()
        D <- as.double(M2[,1,drop=F])
        dim(D) <- c(nrow(M2),1)
        Y <- as.double(M2[,2:4])
        dim(Y) <- c(nrow(M2),3)
        c = single_continuous_fit(D,Y,BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "hill",degree = 4)
        validate_model(c,  "Model: Hill Distribution: Normal" ,  c(6.0730964791569, -5.32930641539168, 39.8852409296944, 3.11573367222951, -0.188967044405645) ,  c(BMD = 24.527770998594, BMDL = 20.2417034249427, BMDU = 29.2708718586858) )
})

test_that("Plots", {
        set.seed(5981)
        M2 <- build_single_continuous_dataset()
        c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "polynomial",degree = 4)
        mle_plot <- plot(c)
        expect_identical(mle_plot$labels$x, "Dose")
        expect_identical(mle_plot$labels$y, "Response")
        expect_identical(mle_plot$labels$title, "Model: Polynomial Distribution: Normal,  Fit Type: MCMC")
        c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                                  distribution = "normal",fit_type="mcmc",model_type = "hill",degree = 4)
        mle_plot <- plot(c)
        expect_identical(mle_plot$labels$x, "Dose")
        expect_identical(mle_plot$labels$y, "Response")
        expect_identical(mle_plot$labels$title, "Model: Hill Distribution: Normal,  Fit Type: MCMC")
})