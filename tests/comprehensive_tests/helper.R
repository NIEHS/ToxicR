library(actuar)

positive_dataset_genes = c("GDF15_33113","A2M_7932","STAC3_9953","ANXA7_8051",
                           "EPHX1_8567","BTG3_33276","APOA1_33150","HSPB1_8847",
                           "ACOT1_7968","MGLL_9227","ACOT2_7969","CCNG1_32302",
                           "ALDH1A1_8022","TM7SF2_10031","ABCC3_7941")

negative_dataset_genes = c("ABCG8_7948",
                           "Ace_34123",
                           "ACTG2_32972",
                           "AK2_33292",
                           "FAM183B_33213",
                           "FSHB_33074",
                           "PROS1_9577",
                           "RT1-A2_33098")

intermediate_dataset_genes = c("ALDH6A1_8027","ATF3_8095","CXCL1_8407",
                               "CYTB_32443","SYPL1_9979","YME1L1_32545")

build_positive_response_dataset <- function(target_gene) {
  mean_data = read.csv("datasets/positive_control_mean.txt", header=FALSE)
  stdev_data = read.csv("datasets/positive_Control_stdev.txt", header=FALSE)
  gene_mean <- unlist(mean_data[grep(target_gene, mean_data[,1]), -1], use.names = FALSE)
  gene_stdev <- unlist(stdev_data[grep(target_gene, stdev_data[,1]), -1], use.names = FALSE)

  response_matrix <- matrix(0, nrow = 9, ncol = 4)
  colnames(response_matrix) <- c("Dose", "Resp", "N", "StDev")
  response_matrix[, 1] <- c(0, 0.156, 0.3125, 0.625, 1.25, 2.5, 5, 10, 20)
  response_matrix[, 2] <- gene_mean
  response_matrix[, 3] <- c(4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0)
  response_matrix[, 4] <- gene_stdev

  return(response_matrix)
}

build_negative_response_dataset <- function(target_gene) {
  mean_data = read.csv("datasets/negative_control_mean.txt", header=FALSE)
  stdev_data = read.csv("datasets/negative_Control_stdev.txt", header=FALSE)
  gene_mean <- unlist(mean_data[grep(target_gene, mean_data[,1]), -1], use.names = FALSE)
  gene_stdev <- unlist(stdev_data[grep(target_gene, stdev_data[,1]), -1], use.names = FALSE)

  response_matrix <- matrix(0, nrow = 9, ncol = 4)
  colnames(response_matrix) <- c("Dose", "Resp", "N", "StDev")
  response_matrix[, 1] <- c(0, 0.156, 0.3125, 0.625, 1.25, 2.5, 5, 10, 20)
  response_matrix[, 2] <- gene_mean
  response_matrix[, 3] <- c(4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0)
  response_matrix[, 4] <- gene_stdev

  return(response_matrix)
}

build_intermediate_response_dataset <- function(target_gene) {
  mean_data = read.csv("datasets/intermediate_control_mean.txt", header=FALSE)
  stdev_data = read.csv("datasets/intermediate_Control_stdev.txt", header=FALSE)
  gene_mean <- unlist(mean_data[grep(target_gene, mean_data[,1]), -1], use.names = FALSE)
  gene_stdev <- unlist(stdev_data[grep(target_gene, stdev_data[,1]), -1], use.names = FALSE)

  response_matrix <- matrix(0, nrow = 9, ncol = 4)
  colnames(response_matrix) <- c("Dose", "Resp", "N", "StDev")
  response_matrix[, 1] <- c(0, 0.156, 0.3125, 0.625, 1.25, 2.5, 5, 10, 20)
  response_matrix[, 2] <- gene_mean
  response_matrix[, 3] <- c(4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0)
  response_matrix[, 4] <- gene_stdev

  return(response_matrix)
}




run_model <- function(M2, model_type, fit_type, distribution){
    model_res <- single_continuous_fit(M2[, 1, drop = F], M2[, 2:4],
                                       BMR_TYPE = "sd", BMR = 1, ewald = T,
                                       distribution = distribution, fit_type = fit_type, model_type = model_type, degree = 4)
    result <- c(model_res$full_model, model_res$parameters, model_res$bmd)
    return(result)
}

model_runner <- function(gene_list, gene_dataset_creator) {
  results <- data.frame(model = character(), distribution = character(), 
                        model_type = character(), fit_type = character(),
                        parameter1 = numeric(),
                        parameter2 = numeric(), parameter3 = numeric(),
                        parameter4 = numeric(), parameter5 = numeric(),
                        parameter6 = numeric(), parameter7 = numeric(),
                        bmd = numeric(), bmdl = numeric(), bmdu = numeric())
  

  for (gene in gene_list){
    M2 <- gene_dataset_creator(gene)
    for (model_type in c("hill", "exp-3", "exp-5", "power", "polynomial", "exp-aerts",
                        "invexp-aerts", "gamma-aerts", "invgamma-aerts", "hill-aerts",
                        "lomax-aerts", "invlomax-aerts", "lognormal-aerts", "logskew-aerts",
                        "invlogskew-aerts", "logistic-aerts", "probit-aerts")){
        for (fit_type in c("laplace", "mle", "mcmc")){
            for (distribution in c("normal", "normal-ncv", "lognormal")){
              if(model_type %in% c("hill", "power", "polynomial") && distribution == "lognormal") {
                # Log-Normal/Hill specification is not presently available.
                # Power-Log-normal models are not allowed. Please choose normal or normal non-constant variance.
                # Polynomial-Log-normal models are not allowed. Please choose normal or normal non-constant variance.
                next
              }
              if(model_type %in% c("exp-aerts",
                        "invexp-aerts", "gamma-aerts", "invgamma-aerts", "hill-aerts",
                        "lomax-aerts", "invlomax-aerts", "lognormal-aerts", "logskew-aerts",
                        "invlogskew-aerts", "logistic-aerts", "probit-aerts") && fit_type == "mle") {
                          # Aerts models are currently not supported with the frequentist approach
                          next
                        }
              result <- single_continuous_fit(M2[, 1, drop = F], M2[, 2:4],
                                       BMR_TYPE = "sd", BMR = 1, ewald = T,
                                       distribution = distribution, fit_type = fit_type, model_type = model_type, degree = 4)
              num_parameters = length(result$parameters)
              bmds <- unlist(result$bmd, use.names=FALSE)
              if (num_parameters == 4) {
                df <- data.frame(model = result$full_model, distribution = distribution, 
                                        model_type = model_type, fit_type = fit_type,
                                        parameter1 = result$parameters[1], 
                                        parameter2 = result$parameters[2], parameter3 = result$parameters[3], 
                                        parameter4 = result$parameters[4], parameter5 = NA, 
                                        parameter6 = NA, parameter7 = NA,
                                        bmd = bmds[1], bmdl = bmds[2], bmdu = bmds[3])
              }
              else if(num_parameters == 5) {
                df <- data.frame(model = result$full_model, distribution = distribution, 
                                        model_type = model_type, fit_type = fit_type,
                                        parameter1 = result$parameters[1], 
                                        parameter2 = result$parameters[2], parameter3 = result$parameters[3], 
                                        parameter4 = result$parameters[4], parameter5 = result$parameters[5], 
                                        parameter6 = NA, parameter7 = NA,
                                        bmd = bmds[1], bmdl = bmds[2], bmdu = bmds[3])
              } else if (num_parameters == 6) {
                df <- data.frame(model = result$full_model, distribution = distribution, 
                                        model_type = model_type, fit_type = fit_type,
                                        parameter1 = result$parameters[1], 
                                        parameter2 = result$parameters[2], parameter3 = result$parameters[3], 
                                        parameter4 = result$parameters[4], parameter5 = result$parameters[5], 
                                        parameter6 = result$parameters[6], parameter7 = NA,
                                        bmd = bmds[1], bmdl = bmds[2], bmdu = bmds[3])
              }
              else if (num_parameters == 7){
                df <- data.frame(model = result$full_model, distribution = distribution, 
                                        model_type = model_type, fit_type = fit_type,
                                        parameter1 = result$parameters[1], 
                                        parameter2 = result$parameters[2], parameter3 = result$parameters[3], 
                                        parameter4 = result$parameters[4], parameter5 = result$parameters[5], 
                                        parameter6 = result$parameters[6], parameter7 = result$parameters[7],
                                        bmd = bmds[1], bmdl = bmds[2], bmdu = bmds[3])
              }
              else{
                print("More than 7 parameters")
                print(num_parameters)
              }
              
              results <- rbind(results, df)
            }
        }
    }

    write.csv(results, file = paste0("positive_control_", gene, "_results.csv"), row.names = FALSE)
    results <- data.frame(model = character(), distribution = character(), 
                        model_type = character(), fit_type = character(),
                        parameter1 = numeric(),
                        parameter2 = numeric(), parameter3 = numeric(),
                        parameter4 = numeric(), parameter5 = numeric(),
                        parameter6 = numeric(), parameter7 = numeric(),
                        bmd = numeric(), bmdl = numeric(), bmdu = numeric())
  
  }
}