context("All models positive control")


set.seed(5981)
gene <- positive_dataset_genes[1]
M2 <- build_positive_response_dataset(gene)
validation_file <- paste0("datasets/", "positive_control_", gene, "_results.csv")

expected_values <- read.csv(validation_file, stringsAsFactors = FALSE)

for (i in seq_len(nrow(expected_values))) {
    test_name <- paste(expected_values[i, "model"], expected_values[i, "distribution"], expected_values[i, "model_type"], expected_values[i, "fit_type"], sep="_")
    test_that(test_name, {
        # extract the expected values for this row
        result <- single_continuous_fit(M2[, 1, drop = F], M2[, 2:4],
                                        BMR_TYPE = "sd", BMR = 1, ewald = T,
                                        distribution = expected_values[i, "distribution"], 
                                        fit_type = expected_values[i, "fit_type"], 
                                        model_type = expected_values[i, "model_type"], 
                                        degree = 4)

        actual_full_name <- result$full_model
        actual_bmd <- result$bmd
        names(actual_bmd) <- NULL

        actual_parameters <- result$parameters
        num_params <- length(actual_parameters) - 1
        idx <- num_params + 5

        expected_full_name <- expected_values[i, "model"]
        expected_bmd <-  unlist(expected_values[i, 12:14], use.names = FALSE)
        expected_parameters <- unlist(expected_values[i, 5:idx], use.names = FALSE)

        expect_equal(actual_full_name, expected_full_name)
        expect_equal(actual_bmd, expected_bmd)
        expect_equal(actual_parameters, expected_parameters)
    })
}

