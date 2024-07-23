context("Single Continuous Model bounds")

test_that("Laplace gamma-efsa bounds are respected", {
    zeros <- c(0.0, 0.0, 0.0, 0.0, 0.0)
    set.seed(1992)
    cont_data <- matrix(0, nrow = 5, ncol = 4)
    colnames(cont_data) <- c("Dose", "Mean", "N", "SD")
    cont_data[, 1] <- c(0, 50, 100, 200, 400)
    cont_data[, 2] <- c(5.26, 5.76, 7.13, 9.24, 9.23)
    cont_data[, 3] <- c(20, 20, 20, 20, 20)
    cont_data[, 4] <- c(2.23, 1.47, 2.47, 2.24, 1.56)
    Y <- cont_data[, 2:4]
    test_new2 <- single_continuous_fit(cont_data[,1],Y, model_type="gamma-efsa", fit_type = 'laplace', distribution='lognormal', alpha = 0.025)
    diag_positive <- diag(test_new2$covariance) > zeros
    expect_equal(diag_positive, c(TRUE, TRUE, TRUE, TRUE, TRUE))
})

test_that("Posterior Probabilities dont differ", {
    cont_data <- matrix(0,nrow=5,ncol=4)
    colnames(cont_data) <- c("Dose","Mean","N","SD")
    cont_data[,1] <- c(0,0.35,1,2.5,5)
    cont_data[,2] <- c(4,4.1,5,7.5,8.2)
    cont_data[,4] <- c(0.58,0.75,1.55,2.44,2.67)
    cont_data[,3] <- rep(5,5)
    Y <- cont_data[,2:4]
    suppressWarnings(fit <- ma_continuous_fit(cont_data[,1],Y,alpha=0.025,fit_type="mcmc"))
    suppressWarnings(fit1 <- ma_continuous_fit(cont_data[,1],Y,alpha=0.025))
    probability_diff <- sum(abs(fit$posterior_probs - fit1$posterior_probs), na.rm = T)
    expect_lte(probability_diff, 0.01)
})

test_that("Negative Hessians don't have posterior probability calculated", {
    cont_data <- matrix(0,nrow=5,ncol=4)
    colnames(cont_data) <- c("Dose","Mean","N","SD")
    cont_data[,1] <- c(0,0.001,0.03,1,20)
    cont_data[,2] <- c(4,8.2,9.1,11.2,12.2)
    cont_data[,4] <- c(0.78,1.46,4.05,3.44,4.67)
    cont_data[,3] <- rep(4,5)
    Y <- cont_data[,2:4]
    fit_types <- c("mcmc", "mle", "laplace")
    for (fit_type in fit_types) {
        suppressWarnings(fit <- ma_continuous_fit(cont_data[,1],Y,alpha=0.025,fit_type=fit_type))
        all_model_names <- names(fit$posterior_probs)
        na_model_names <- names(which(is.na(fit$posterior_probs)))
        model_name_prefix <- "Indiv_"
        for (name in na_model_names) {
            new_name <- paste0(model_name_prefix, name)
            determinant <- det(fit[[new_name]]$covariance)
            is_bmd_na <- is.na(fit[[new_name]]$bmd[1])
            expect_equal(TRUE, is_bmd_na || determinant < 0.0)
        }
    }
})