library(actuar)

build_ma_dataset_2 <- function() {
        Hill.p <- rbind(
                c(481, -250.3, 70, 3.3),
                c(481, -250.3, 40, 1.3),
                c(481, -250.2, 15, 1.1),
                c(481, -250.3, 50, 4),
                c(10.58, 9.7, 70, 3.5),
                c(10.58, 9.7, 25, 3),
                c(10.58, 9.7, 15, 2),
                c(10.58, 9.7, 50, 4)
        )
        hill <- data.frame(
                a = Hill.p[, 1], b = Hill.p[, 2],
                c = Hill.p[, 3], d = Hill.p[, 4]
        )


        doses <- rep(c(0, 6.25, 12.5, 25, 50, 80, 100, 150), each = 10)

        mean <- ToxicR:::.cont_hill_f(as.numeric(hill[1, ]), doses)
        set.seed(2020)
        y <- rinvgauss(length(mean), mean, 45528.14)
        return(list(doses = doses, y = y))
}

build_single_continuous_dataset <- function() {
        M2 <- matrix(0, nrow = 5, ncol = 4)
        colnames(M2) <- c("Dose", "Resp", "N", "StDev")
        M2[, 1] <- c(0, 25, 50, 100, 200)
        M2[, 3] <- c(20, 20, 19, 20, 20)
        M2[, 2] <- c(6, 5.2, 2.4, 1.1, 0.75)
        M2[, 4] <- c(1.2, 1.1, 0.81, 0.74, 0.66)
        M2
}

build_single_dichotomous_dataset_2 <- function() {
        mData <- matrix(c(
                0, 39, 297,
                0.00098, 24, 90,
                0.0098, 32, 87,
                0.098, 136, 148
        ),
        nrow = 4, ncol = 3, byrow = T
        )
        return(mData)
}


validate_model2 <- function(model, name, parameters, bmd_estimates,
                            gof) {
        expect_equal(name, model$full_model)
        expect_equal(parameters, model$parameters, tolerance = 10e-3)
        expect_equal(setNames(bmd_estimates, c("BMD", "BMDL", "BMDU")), model$bmd, tolerance = 10e-3)
        A <- summary(model)
        expect_equal(as.numeric(A$GOF), gof, tolerance = 10e-3)
}

build_single_dichotomous_dataset <- function() {
        mData <- matrix(c(
                0, 2, 50,
                1, 2, 50,
                3, 10, 50,
                16, 18, 50,
                32, 18, 50,
                33, 17, 50
        ), nrow = 6, ncol = 3, byrow = T)
        mData
}

build_single_dichotomous_dataset2 <- function(){
  mData <- matrix(c(0, 2,50,
                    1, 2,50,
                    3, 10, 50,
                    16, 13,50,
                    32, 18,50,
                    33, 20,50),nrow=6,ncol=3,byrow=T)
  mData
}

validate_model <- function(model, name, parameters, bmd_estimates, tolerance = 10e-2) {
        expect_equal(name, model$full_model)
        expect_equal(parameters, model$parameters, tolerance = tolerance)
        expect_equal(setNames(bmd_estimates, c("BMD", "BMDL", "BMDU")), model$bmd, tolerance = tolerance)
        # show(model$full_model)
        # show(model$parameters)
        # show(model$bmd)
}

generate_validation_code_single <- function(c) {
        cat("\n")
        if ("fitted_model" %in% names(c)) {
                cat("validate_model(c, ", paste0("\"", c$fitted_model$full_model, "\""), ", ", paste(list(c$fitted_model$parameters), sep = ", "), ", ", paste(list(c$bmd), sep = ", "), ")\n")
        } else {
                cat("validate_model(c, ", paste0("\"", c$full_model, "\""), ", ", paste(list(c$parameters), sep = ", "), ", ", paste(list(c$bmd), sep = ", "), ")\n")
        }
}


generate_validation_code <- function(AA) {
        cat("\n")
        for (i in 1:(length(AA) - 3)) {
                cat("validate_model(", paste0("AA$Individual_Model_", i), ", ", paste0("\"", AA[[i]]$full_model, "\""), ", ", paste(list(AA[[i]]$parameters), sep = ", "), ", ", paste(list(AA[[i]]$bmd), sep = ", "), ")\n")
        }
}

build_ma_dataset <- function() {
        Hill.p <- rbind(
                c(481, -250.3, 70, 3.3),
                c(481, -250.3, 40, 1.3),
                c(481, -250.2, 15, 1.1),
                c(481, -250.3, 50, 4),
                c(10.58, 9.7, 70, 3.5),
                c(10.58, 9.7, 25, 3),
                c(10.58, 9.7, 15, 2),
                c(10.58, 9.7, 50, 4)
        )
        hill <- data.frame(
                a = Hill.p[, 1], b = Hill.p[, 2],
                c = Hill.p[, 3], d = Hill.p[, 4]
        )


        doses <- rep(c(0, 6.25, 12.5, 25, 50, 100), each = 10)
        dosesq <- rep(c(0, 6.25, 12.5, 25, 50, 100), each = 30)

        mean <- ToxicR:::.cont_hill_f(as.numeric(hill[6, ]), doses)
        y <- rinvgauss(length(mean), mean, 18528.14)
        return(list(doses = doses, y = y))
}

build_model_list <- function(y) {
        model_listA <- data.frame(
                model_list = c(rep("hill", 2), rep("exp-3", 3), rep("exp-5", 3), rep("power", 2)),
                distribution_list = c(
                        "normal", "normal-ncv", rep(c("normal", "normal-ncv", "lognormal"), 2), "normal",
                        "normal-ncv"
                )
        )
        model_list <- list()
        for (i in 1:nrow(model_listA)) {
                t_prior <- ToxicR:::.bayesian_prior_continuous_default(model_listA$model_list[i], model_listA$distribution_list[i])
                if (model_listA$distribution_list[i] == "lognormal") {
                        t_prior$priors[nrow(t_prior$priors), 2] <- log(var(log(y)))
                } else {
                        if (model_listA$distribution_list[i] == "normal") {
                                t_prior$priors[nrow(t_prior$priors), 2] <- log(var(y))
                        } else {
                                t_prior$priors[nrow(t_prior$priors), 2] <- log(mean(y) / var(y))
                        }
                }

                model_list[[i]] <- create_continuous_prior(t_prior, model_listA$model_list[i], model_listA$distribution_list[i])
        }
        model_list
}
