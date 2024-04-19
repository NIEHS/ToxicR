# Copyright 2021  NIEHS <matt.wheeler@nih.gov>
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#' 733 unique dose-response datasets
#'
#' A dataset containing 733 dichotomous dose-response studies that were involved in
#' regulatory risk assessment.
#'
#' @format A data frame with 2727 rows and 11 variables:
#' \describe{
#'   \item{ID}{-The study ID in the database.}
#'   \item{chemical}{-Name of the Chemical in the study.}
#'   \item{data.source}{-Source of the risk assessment data.}
#'   \item{CASRN}{-Chemical's CASRN}
#'   \item{dose}{-Dose spacing of the study using the original study.}
#'   \item{r.dose}{-Doses of the experiment relative to 1 being the maximum dose tested.}
#'   \item{n}{-Number of animals on test.}
#'   \item{obs}{-Number of adverse events.}
#'   \item{organ}{-Organ impacted.}
#'   \item{effect}{-Type of adverse effect.}
#'   \item{study.source}{-Publication related to the experiment.}
#' }
#' More information at: \doi{10.1111/risa.13218}
"dichotomousDR"

#' Short term terminal body-weight data from NTP Report 599
#'
#' This dataset contains terminal body-weight data for male and
#' female rats for the technical report TR-599: Sodium Tungstate Dihydrate.
#'
#' @format A data frame with 120 rows and 4 variables:
#' \describe{
#'   \item{Dose_Group}{-The dose group for the observation.}
#'   \item{dose}{-The dose in mg/L }
#'   \item{sex}{-Animal's Sex}
#'   \item{weight}{-Terminal body-weight}
#' }
#' For more information see: \doi{10.22427/NTP-DATA-TR-599}
"ntp_weight_data"

#' Long term Thyroid Adenoma data from NTP Report 599
#'
#' This dataset contains Thyroid Adenoma data  for
#' female rats for the technical report TR-599: Sodium Tungstate Dihydrate.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{treatment}{-The dose group for the observation.}
#'   \item{days_on_study}{-Number of days on the study 730 is the max.}
#'   \item{adenoma}{- Thyroid Adenoma (Yes/No) (1/0).}
#'   \item{dose}{-The dose in mg/L}
#' }
#' For more information see: \doi{10.22427/NTP-DATA-TR-599}
"ntp_599_female"

#' Clinical Chemistry data from NTP Report 599
#'
#' This dataset contains clinical chemistry data  for
#' all rats in the short term 90-day study.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{concentration}{-The dose group for the observation.}
#'   \item{sex}{- Male/Female.}
#'   \item{response}{- Response variable}
#'   \item{response_type}{- The type of response measured}
#' }
#' For more information see: \doi{10.22427/NTP-DATA-TR-599}
"ntp_599_hemotology"

## ----------------------
## 	POLYK-TEST
## ----------------------
#' @title Poly-k test
#' This function implements the NTP's polyK trend test.
#' @param dose An equation of the form \eqn{Y \sim X.} Here the variable
#' \eqn{Y} is the response of interest, and \eqn{X} represents discrete experimental
#' conditions. For example, if weight is the dependent variable, and you are
#' interested in looking at the trend across sex one would have 'weight ~ sex'.
#' @param tumor A data frame with column names in the formula.
#' @param daysOnStudy The name of the variable containing the doses in the data frame \eqn{data}.
#' It is expected multiple doses for each of the experimental conditions \eqn{X}.
#' @return The results of a Williams trend test for each level in dose_name.
#' More information on this procedure at: \doi{10.2307/2531856} and \doi{10.2307/2532200}
#' This procedure returns a vector of three p-values for the poly-1.5, poly-3, and poly-6 test
#' respectively.
#' @examples
#' ntp_polyk(ntp_599_female$dose, ntp_599_female$adenoma, ntp_599_female$days_on_study)
ntp_polyk <- function(dose, tumor, daysOnStudy) {
  if (sum(tumor > 1) > 0) {
    stop("Tumors need to be a 0 or 1")
  }
  if (sum(tumor < 0) > 0) {
    stop("Tumors need to be a 0 or 1")
  }
  if (sum(daysOnStudy < 0) > 0) {
    stop("Can not have negative days on study.")
  }

  if ((length(dose) != length(tumor)) ||
    (length(tumor) != length(daysOnStudy))) {
    stop("All arrays are not of the same length.")
  }
  if ((sum(is.na(dose)) + sum(is.na(tumor)) + sum(is.na(daysOnStudy))) > 0) {
    stop("There is an NA in the data.")
  }
  result <- .polykCPP(dose, tumor, daysOnStudy)

  result <- as.matrix(result)
  row.names(result) <- c("Poly 1.5", "Poly-3", "Poly-6")
  class(result) <- "ntp.polyk"
  return(result)
}

.print_polyk_ntp <- function(x, ...) {
  result <- x
  cat("The results of the Poly-K test for trend.\n")
  cat(sprintf("Poly-1.5 P-value = %1.4f\n", result[1]))
  cat(sprintf("Poly-3   P-value = %1.4f\n", result[2]))
  cat(sprintf("Poly-6   P-value = %1.4f\n", result[3]))
}
## -----------------------------------------------------------
## JONCKHEERE'S TEST
## ----------------Changelog----------------------------------
## Released:
## [1.1.0] - 08/07/2019. Jira ticket:CEBSPROV-5301
## Changed
## Fixes for Kendall test. Set Exact = FALSE
## Changed to make it based upon a general formula specified in
## formula. To do this, we assume that data is a data frame.
## As a default, "dose_name", is set to the column header "dose"
## -----------------------------------------------------------
#' @title ntp_jonckeere
#' Jonckherre's test for significant differences from background dose
#' @param formula An equation of the form \eqn{Y \sim X.} Here the variable
#' \eqn{Y} is the response of interest, and \eqn{X} represents discrete experimental
#' conditions. For example, if weight is the dependent variable, and you are
#' interested in looking at the trend across sex one would have 'weight ~ sex'.
#' @param data A data frame with column names in the formula.
#' @param dose_name The name of the variable containing the doses in the data frame \eqn{data}.
#' It is expected multiple doses for each of the experimental conditions \eqn{X}.
#' @param pair  The type of test used for pairwise comparison. It can either be
#' "Williams" or "Shirley"
#' @return The results of a global test for difference from background.
#' @examples
#'
#' ntp_jonckeere(
#'      response ~ sex + response_type, 
#'      data = ntp_599_hemotology, 
#'      dose_name = "concentration"
#' )
ntp_jonckeere <- function(formula, data, dose_name = "dose", pair = "Williams") {
  if (!is.data.frame(data)) {
    stop("The data do not appear to be in the data frame format.")
  }
  if (sum(colnames(data) %in% dose_name) == 0) {
    stop(sprintf("The specified dose column '%s' does not exist in the data frame.", dose_name))
  }
  if (!(pair %in% c("Williams", "Shirley"))) {
    stop("Variable 'pair' must be either 'Williams' or 'Shirley'")
  }


  jonck <- NULL
  ## KRS - added "numeric_value" on the left hand side
  jonck_list <- as.data.frame(summaryBy(formula, data = data, FUN = length))
  jonck_list[is.na(jonck_list)] <- ""

  ## may create WARNINGS when ties are present
  analysis_var <- as.character(unlist(formula[[2]]))
  temp_colnames <- unlist(c(unlist(colnames(jonck_list)), dose_name, as.character(unlist(formula[[2]]))))
  temp <- colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[, temp == TRUE])

  for (i in 1:nrow(jonck_list))
  {
    temp_names <- rep(NA, ncol(jonck_list))
    for (j in 1:length(temp_names)) {
      temp_names[j] <- unlist(jonck_list[i, j])
    }

    ## KRS - changed "phase" to "phasetype"
    temp_dd <- temp_d
    temp_dd[is.na(temp_dd)] <- ""

    for (j in 1:length(jonck_list)) {
      myt <- which(colnames(temp_d) == colnames(jonck_list)[j])

      if (length(myt) > 0) {
        temp_dd <- temp_dd[as.character(temp_dd[, myt]) == as.character(temp_names[j]), ]
      }
    }
    tempdata <- temp_dd
    tidx <- which(colnames(tempdata) == dose_name)
    # print(tempdata)
    ## make sure there is more than ONE record for gender i
    if (nrow(tempdata) > 1) {
      ## make sure variance is not zero (creates error) and control is present
      if ((length(unique(tempdata[, analysis_var])) > 1) & (0 %in% as.numeric(tempdata[, tidx]))) {
        stat <- cor.test(as.numeric(tempdata[, tidx]), tempdata[, analysis_var], method = "kendall", exact = FALSE)
        jline <- c(temp_names, stat$estimate, stat$p.value)
        jonck <- rbind(jonck, jline)
      }
    }
  }

  ## check for existence
  if (!is.null(jonck)) {
    rownames(jonck) <- NULL
    jonck <- as.data.frame(jonck)
    names(jonck) <- c(colnames(jonck_list), "tau", "pvalue")
    jonck$tau <- as.numeric(as.character(jonck$tau))
    jonck$pvalue <- as.numeric(as.character(jonck$pvalue))

    ## need to remove nulls
    jonck[is.na(jonck)] <- ""

    temp_col <- sprintf("%s.length", as.character(formula[[2]]))
    temp_idx <- which(temp_col == colnames(jonck))
    jonck <- jonck[, -temp_idx]

    if (pair == "Williams") {
      jonck$mult_comp_test <- ifelse(jonck$pvalue <= .01, "WILLIAMS", "DUNNETT")
    }
    if (pair == "Shirley") {
      jonck$mult_comp_test <- ifelse(jonck$pvalue <= .01, "SHIRLEY", "DUNN")
    }
  }
  return(jonck)
}

.compute_crit_williams <- function(william_test_data, dose_name, formulaV) {
  dof <- NULL

  t_idx <- which(colnames(william_test_data) == dose_name)
  william_test_data[, t_idx] <- as.numeric(william_test_data[, t_idx])
  william_test_data$crit05 <- NA
  william_test_data$crit01 <- NA

  for (k in 1:nrow(william_test_data)) {
    ## CONTROL GROUP
    if (william_test_data[k, t_idx] == 0) ## should this be datatemp or datatemp?
      {
        rm
        william_test_data$crit05 <- FALSE
        william_test_data$crit01 <- FALSE
      } else if (william_test_data[k, t_idx] != 0 & (william_test_data$dof[k] %in% .will005$dof)) {
      col1 <- paste("w1crit", k, sep = "")
      col5 <- paste("w5crit", k, sep = "")
      adj1 <- paste("w1adj", k, sep = "")
      adj5 <- paste("w5adj", k, sep = "")

      w1crit <- subset(.will005, dof == william_test_data$dof[k])[, c(col1)]
      w1adj <- subset(.will005, dof == william_test_data$dof[k])[, c(adj1)]
      w5crit <- subset(.will025, dof == william_test_data$dof[k])[, c(col5)]
      w5adj <- subset(.will025, dof == william_test_data$dof[k])[, c(adj5)]

      dofactor <- ((william_test_data$dof[k] - lowdof) / (highdof - lowdof))
      temp_name <- sprintf("%s.length", formulaV)
      t_idx2 <- which(colnames(william_test_data) == temp_name)

      con_num <- william_test_data[william_test_data[, t_idx] == 0, ][, t_idx2]
      trt_num <- william_test_data[k, t_idx2]


      william_test_data$crit01[k] <- w1crit - (.1 * w1adj * (1 - (trt_num / con_num)))
      william_test_data$crit05[k] <- w5crit - (.1 * w5adj * (1 - (trt_num / con_num)))
    } else if (william_test_data[k, t_idx] != 0 & !(william_test_data$dof[k] %in% .will005$dof)) {
      col1 <- paste("w1crit", k, sep = "")
      col5 <- paste("w5crit", k, sep = "")
      adj1 <- paste("w1adj", k, sep = "")
      adj5 <- paste("w5adj", k, sep = "")

      ## get lower bound from table
      lowdof <- max(.will005$dof[william_test_data$dof[k] > .will005$dof])

      low.w1crit <- subset(.will005, dof == lowdof)[, c(col1)]
      low.w1adj <- subset(.will005, dof == lowdof)[, c(adj1)]
      low.w5crit <- subset(.will025, dof == lowdof)[, c(col5)]
      low.w5adj <- subset(.will025, dof == lowdof)[, c(adj5)]

      ## get upper bound from table
      highdof <- min(.will005$dof[william_test_data$dof[k] < .will005$dof])

      high.w1crit <- subset(.will005, dof == highdof)[, c(col1)]
      high.w1adj <- subset(.will005, dof == highdof)[, c(adj1)]
      high.w5crit <- subset(.will025, dof == highdof)[, c(col5)]
      high.w5adj <- subset(.will025, dof == highdof)[, c(adj5)]

      dofactor <- ((william_test_data$dof[k] - lowdof) / (highdof - lowdof))
      temp_name <- sprintf("%s.length", formulaV)
      t_idx2 <- which(colnames(william_test_data) == temp_name)

      con_num <- william_test_data[william_test_data[, t_idx] == 0, ][, t_idx2]
      trt_num <- william_test_data[k, t_idx2]

      william_test_data$crit01[k] <- (low.w1crit - (dofactor * (low.w1crit - high.w1crit))) - (.01 * low.w1adj * (1 - (trt_num / con_num)))
      william_test_data$crit05[k] <- (low.w5crit - (dofactor * (low.w5crit - high.w5crit))) - (.01 * low.w5adj * (1 - (trt_num / con_num)))
    }
  }
  return(william_test_data)
}

## ----------------------
## 	WILLIAM'S TEST
## ----------------------
#' Williams Trend test for
#' @title Wiliam's trend test
#' @param formula An equation of the form \eqn{Y \sim X.} Here the variable
#' \eqn{Y} is the response of interest, and \eqn{X} represents discrete experimental
#' conditions. For example, if weight is the dependent variable, and you are
#' interested in looking at the trend across sex one would have 'weight ~ sex'.
#' @param data A data frame with column names in the formula.
#' @param dose_name The name of the variable containing the doses in the data frame \eqn{data}.
#' It is expected multiple doses for each of the experimental conditions \eqn{X}.
#' @return The results of a Williams trend test for each level in \eqn{dose_name}.
#' For more information on the Williams trend test: \doi{10.2307/2528930}
#' #' \itemize{
#'  \item \code{X}: this represents all the class objects on the right hand side of \eqn{ Y \sim X} above.
#'  \item \code{dose:} the dose groups relative to control.
#'  \item \code{willStat}: Value of the Shirley test statistic.
#'  \item \code{mult_comp_signif}: Test's significance as 0, 1, or 2 which is not-significant,
#'  significant at the 0.05% level and significant at the 0.01% level.
#'  \item \code{mult_comp_test}: The type of test, i.e. "William"
#' }
#' @examples
#'
#' a <- ntp_williams(weight ~ sex, data = ntp_weight_data)
#' summary(a)
ntp_williams <- function(formula, data, dose_name = "dose") {
  mult_comp_test <- NULL
  data[, c(dose_name)] <- as.numeric(unlist(data[, c(dose_name)]))

  temp_str <- strsplit(as.character(formula)[3], " ")[[1]]
  temp_str <- temp_str[temp_str != "+"]
  data <- data[order(data[, c(dose_name)]), ]

  for (ii in 1:length(temp_str)) {
    data <- data[order(data[, temp_str[ii]]), ]
  }

  if (!(dose_name %in% colnames(data))) {
    stop(sprintf("Dose name %s does not appear in the data frame.", dose_name))
  }

  jonck_data <- ntp_jonckeere(formula,
    dose_name = dose_name,
    data = data, pair = "Williams"
  )

  ## loop through all groups flagged as WILLIAM in jonck
  william <- subset(jonck_data, mult_comp_test == "WILLIAMS")
  will_results <- NULL
  will_results2 <- NULL

  temp_resp_name <- unlist(formula[[2]])

  temp_colnames <- unlist(c(unlist(colnames(william)), dose_name, as.character(unlist(formula[[2]]))))
  temp <- colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[, temp == TRUE])

  if (nrow(william) > 0) {
    for (w in 1:nrow(william)) {
      temp_names <- rep(NA, ncol(william))
      for (j in 1:length(temp_names)) {
        temp_names[j] <- unlist(william[w, j])
      }
      ## KRS - changed "phase" to "phasetype"
      temp_dd <- temp_d
      temp_dd[is.na(temp_dd)] <- ""
      # subset the data based upon the columns that exist in the formula
      for (j in 1:length(william)) {
        myt <- which(colnames(temp_d) == colnames(william)[j])

        if (length(myt) > 0) {
          temp_dd <- temp_dd[as.character(temp_dd[, myt]) == as.character(temp_names[j]), ]
        }
      }

      datatemp <- temp_dd
      temp_idx <- which(colnames(datatemp) == dose_name)
      datatemp <- datatemp[order(datatemp[, temp_idx]), ]
      ## get william-ized dose means
      ## direction of smoothing dependent on Jonckheere output
      formula_temp <- sprintf("%s ~ %s + %s", as.character(formula[[2]]), as.character(dose_name), as.character(formula)[3])
      wmeans <- summaryBy(as.formula(formula_temp), data = datatemp, FUN = c(mean, length))
      temp_idx2 <- which(colnames(wmeans) == dose_name)
      wmeans <- wmeans[order(wmeans[, temp_idx2]), ]
      temp_value <- paste(temp_resp_name, ".mean", sep = "")
      smeans <- wmeans[, temp_value]
      temp_value <- paste(temp_resp_name, ".length", sep = "")
      slength <- wmeans[, temp_value]
      # smeans  <- wmeans$numeric_value.mean
      # slength <- wmeans$numeric_value.length

      smeans <- smeans[-1] ## remove control group
      slength <- slength[-1] ## remove control group

      trend <- ifelse(william$tau[w] < 0, william$pvalue[w] * -1, william$pvalue[w])
      ## set comparison direction based on JONCK trend result
      direction <- ifelse(trend < 0, "decreasing", "increasing")

      ## need more than 1 trt grp for smoothing
      if (length(smeans) > 1) {
        if (direction == "decreasing") {
          for (i in 1:(length(smeans) - 1))
          {
            if (smeans[i] < smeans[i + 1]) ## smoothing required
              {
                if (i == 1) ## FIRST ITEM
                  {
                    tempSmean <- (smeans[i] * (slength[i] / (slength[i] + slength[i + 1]))) + (smeans[i + 1]) * (slength[i + 1] / (slength[i] + slength[i + 1]))
                    smeans[i] <- tempSmean
                    smeans[i + 1] <- tempSmean
                  } else

                if (i > 1) ## NOT FIRST ITEM
                  {
                    tempSmean <- (smeans[i] * (slength[i] / (slength[i] + slength[i + 1]))) + (smeans[i + 1]) * (slength[i + 1] / (slength[i] + slength[i + 1]))

                    s_count <- 0 ## initialize counter for previous consecutive smoothed means
                    for (j in ((i - 1):1)) ## loop back to find number of elements in smoothing
                    {
                      if (tempSmean > smeans[j]) {
                        s_count <- s_count + 1 ## count once for each consecutive previous "smooth"
                        if (s_count > 0) {
                          tempSmean <- 0
                          for (k in (i - s_count:i + 1)) ## recalculate tempSmean
                          {
                            tempSmean <- tempSmean + smeans[k] * (slength[k] / sum(slength[i - s_count:i + 1]))
                          }
                        }
                      } else {
                        break
                      }
                    }

                    denom <- 0 ## this is denominator for weighting by sample size
                    for (k in (i + 1):(i - s_count)) ## get count from previous smoothing
                    {
                      denom <- denom + slength[k]
                    }


                    tempSmean <- 0
                    for (k in (i + 1):(i - s_count)) ## go back through and get smoothed mean(s)
                    {
                      tempSmean <- tempSmean + smeans[k] * slength[k] / denom
                    }


                    smeans[i] <- tempSmean
                    smeans[i + 1] <- tempSmean
                    if (s_count > 0) {
                      for (k in 1:s_count)
                      {
                        smeans[i - k] <- tempSmean
                      }
                    }
                  }
              }
          }
        } else ## DIRECTION BREAK
        {
          for (i in 1:(length(smeans) - 1))
          {
            if (smeans[i] > smeans[i + 1]) ## smoothing required
              {
                if (i == 1) ## FIRST ITEM
                  {
                    tempSmean <- (smeans[i] * (slength[i] / (slength[i] + slength[i + 1]))) + (smeans[i + 1]) * (slength[i + 1] / (slength[i] + slength[i + 1]))
                    smeans[i] <- tempSmean
                    smeans[i + 1] <- tempSmean
                  } else

                if (i > 1) ## NOT FIRST ITEM
                  {
                    tempSmean <- (smeans[i] * (slength[i] / (slength[i] + slength[i + 1]))) + (smeans[i + 1]) * (slength[i + 1] / (slength[i] + slength[i + 1]))

                    s_count <- 0 ## initialize counter for previous consecutive smoothed means
                    for (j in ((i - 1):1)) ## loop back to find number of elements in smoothing
                    {
                      if (tempSmean < smeans[j]) {
                        s_count <- s_count + 1 ## count once for each consecutive previous "smooth"
                        if (s_count > 0) {
                          tempSmean <- 0
                          for (k in (i - s_count:i + 1)) ## recalculate tempSmean
                          {
                            tempSmean <- tempSmean + smeans[k] * (slength[k] / sum(slength[i - s_count:i + 1]))
                          }
                        }
                      } else {
                        break
                      }
                    }

                    denom <- 0 ## this is denominator for weighting by sample size
                    for (k in (i + 1):(i - s_count)) ## get count from previous smoothing
                    {
                      denom <- denom + slength[k]
                    }


                    tempSmean <- 0
                    for (k in (i + 1):(i - s_count)) ## go back through and get smoothed mean(s)
                    {
                      tempSmean <- tempSmean + smeans[k] * slength[k] / denom
                    }


                    smeans[i] <- tempSmean
                    smeans[i + 1] <- tempSmean
                    if (s_count > 0) {
                      for (k in 1:s_count)
                      {
                        smeans[i - k] <- tempSmean
                      }
                    }
                  }
              }
          }
        }
      }
      temp_value <- paste(temp_resp_name, ".mean", sep = "")
      smeans <- c(wmeans[1, temp_value], smeans)

      # smeans <- c(wmeans$numeric_value.mean[1], smeans)	## pre-pend control mean to beginning of smoothed means
      wmeans <- cbind(wmeans, smeans) ## combine with means info

      ##
      temp_value <- sprintf("%s.length", as.character(formula[[2]]))
      temp_idx4 <- which(colnames(wmeans) == temp_value)
      ## simplify dof calcs ... if errors try old method above
      dof1 <- sum(wmeans[, 5])
      dof2 <- nrow(wmeans)
      wmeans$dof <- (dof1 - dof2)
      formula_text <- sprintf("%s ~ %s", as.character(formula[[2]]), dose_name)
      se <- summaryBy(as.formula(formula_text), data = datatemp, FUN = c(sd, length))
      temp_value <- sprintf("%s.sd", as.character(formula[[2]]))
      temp_idx4 <- which(colnames(se) == temp_value)
      temp_value <- sprintf("%s.length", as.character(formula[[2]]))
      temp_idx5 <- which(colnames(se) == temp_value)
      se$sterr <- se[, temp_idx4] / (se[, temp_idx5]^.5)
      temp_idx4 <- which(colnames(se) == dose_name)
      temp_idx5 <- which(colnames(datatemp) == dose_name)
      mse <- 0
      for (i in 1:nrow(se)) {
        tempS <- sum(datatemp[, temp_idx5] == se[i, temp_idx4])
        temp <- se$sterr[i]^2 * tempS * (tempS - 1)
        mse <- mse + temp
      }

      mse <- mse / (nrow(datatemp) - length(unique(datatemp[, temp_idx5])))

      ## create WILLIAMS TEST STATISTIC
      willStat <- NULL

      temp_value <- sprintf("%s.length", as.character(formula[[2]]))
      temp_idx4 <- which(colnames(wmeans) == temp_value)

      for (i in 2:nrow(wmeans)) {
        control_num <- wmeans[1, temp_idx4]
        control_mean <- wmeans$smeans[1]
        test_num <- wmeans[i, temp_idx4]
        test_mean <- wmeans$smeans[i]

        willstat <- ifelse(direction == "decreasing", (control_mean - test_mean) / ((mse * ((1 / test_num) + (1 / control_num)))^.5), (test_mean - control_mean) / ((mse * ((1 / test_num) + (1 / control_num)))^.5))
        willStat <- c(willStat, willstat)
      }

      willStat <- c("control", willStat)
      wmeans$willStat <- willStat
      wmeans <- .compute_crit_williams(wmeans, dose_name = dose_name, formulaV = as.character(formula[[2]])) # compute the critical values for the Williams Trend Test
      will_results <- rbind(will_results, wmeans)
    }

    twill <- which(!(colnames(william) %in% c("tau", "pvalue", "mult_comp_test")))
    idx <- c(which(colnames(will_results) %in% colnames(william)[twill]), which(colnames(will_results) == "dose"))
    # way too complicated loop to do something simple *sigh*
    for (ii in length(idx):1) {
      reorder <- sort(will_results[, idx[ii]], index = TRUE)$ix
      will_results <- will_results[reorder, ]
    }

    ## ----------------------------------------------------------------------
    ## convert williams statistic into p-value based on SAS crit levels
    ## read in critical tables
    ## ----------------------------------------------------------------------
    ## remove control
    will_results2 <- subset(will_results, willStat != "control")

    will_results2$willStat <- as.numeric(as.character(will_results2$willStat))
    will_results2$crit05 <- as.numeric(as.character(will_results2$crit05))
    will_results2$crit01 <- as.numeric(as.character(will_results2$crit01))

    ## determine how many asterisks each row deserves
    will_results2$mult_comp_signif <- ifelse(will_results2$willStat >= will_results2$crit01, 2,
      ifelse(will_results2$willStat >= will_results2$crit05, 1, 0)
    )

    will_results2$mult_comp_test <- "WILLIAMS"
    ta1 <- sprintf("%s.mean", as.character(formula[[2]]))
    ta2 <- sprintf("%s.length", as.character(formula[[2]]))
    t_idx <- which(!(colnames(will_results2) %in% c(ta1, ta2, "crit05", "crit01", "smeans", "dof")))
    will_results2 <- will_results2[, t_idx]
  }

  class(will_results2) <- "ntp.williams"
  return(will_results2)
}

## ------------------------
## DUNN'S TEST
## ------------------------
#' @title ntp_dunn Dunn's test
#' @param formula An equation of the form \eqn{Y \sim X.} Here the variable
#' \eqn{Y} is the response of interest, and \eqn{X} represents discrete experimental
#' conditions. For example, if weight is the dependent variable, and you are
#' interested in looking at the trend across sex one would have 'weight ~ sex'.
#' @param data A data frame with column names in the formula.
#' @param dose_name The name of the variable containing the doses in the data frame \eqn{data}.
#' It is expected multiple doses for each of the experimental conditions \eqn{X}.
#' @return The results of a Dunn's  test for each level in \eqn{dose_name}.
#' @examples
#'
#' a <- ntp_dunn(response ~ sex + response_type,
#'   data = ntp_599_hemotology,
#'   dose_name = "concentration"
#' )
#' summary(a)
ntp_dunn <- function(formula, data, dose_name = "dose") {
  mult_comp_test <- numTies <- NULL
  data[, c(dose_name)] <- as.numeric(unlist(data[, c(dose_name)]))
  temp_str <- strsplit(as.character(formula)[3], " ")[[1]]
  temp_str <- temp_str[temp_str != "+"]
  data <- data[order(data[, c(dose_name)]), ]

  for (ii in 1:length(temp_str)) {
    data <- data[order(data[, temp_str[ii]]), ]
  }

  if (!(dose_name %in% colnames(data))) {
    stop(sprintf("Dose name %s does not appear in the data frame.", dose_name))
  }

  # first do Jonckheere's Test and subset all the 'DUNN' tests
  jonck_data <- ntp_jonckeere(formula,
    dose_name = dose_name,
    data = data, pair = "Shirley"
  )
  dunn <- subset(jonck_data, mult_comp_test == "DUNN")
  dunn_results <- NULL


  temp_colnames <- unlist(c(unlist(colnames(dunn)), dose_name, as.character(unlist(formula[[2]]))))
  temp <- colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[, temp == TRUE])
  if (nrow(dunn) == 0) {
    warning("The Jonckherre test did not suggest the test be performed. Returning a NULL dataset.")
  }
  ## loop through all groups flagged as DUNN in jonck
  if (nrow(dunn) > 0) {
    for (d in 1:nrow(dunn)) {
      temp_names <- rep(NA, ncol(dunn))
      for (j in 1:length(temp_names)) {
        temp_names[j] <- unlist(dunn[d, j])
      }

      ## KRS - changed "phase" to "phasetype"
      temp_dd <- temp_d
      temp_dd[is.na(temp_dd)] <- ""
      # subset the data based upon the columns that exist in the formula
      for (j in 1:length(dunn)) {
        myt <- which(colnames(temp_d) == colnames(dunn)[j])

        if (length(myt) > 0) {
          temp_dd <- temp_dd[as.character(temp_dd[, myt]) == as.character(temp_names[j]), ]
        }
      }
      ex <- temp_dd
      dose_idx <- which(colnames(ex) == dose_name)
      nval_idx <- which(colnames(ex) == as.character(unlist(formula[[2]])))
      ex[, dose_idx] <- as.numeric(ex[, dose_idx])

      ## make sure there are multiple doses to compare, and control is present
      if ((length(unique(ex[, dose_idx])) > 1) & (0 %in% unique(ex[, dose_idx]))) {

        ## rank among each sex/test_name/phase/etc.
        ties <- as.data.frame(table(table(ex[, nval_idx])))

        if (nrow(ties) > 0) {
          names(ties) <- c("numTies", "count")
          ties$numTies <- as.integer(as.character(ties$numTies))

          ## remove singletons
          ties <- subset(ties, numTies != 1)

          ## tie adjustment correction
          correction <- 0
          if (nrow(ties) > 0) {
            for (i in 1:nrow(ties))
            {
              correction <- correction + (ties$count[i] * (ties$numTies[i]^3 - ties$numTies[i]))
            }
          }
        }
        ## rankings ... lowest value gets rank of 1
        ranks <- as.data.frame(table(ex[, nval_idx]))
        names(ranks) <- c("value", "count")
        ranks$value <- as.numeric(as.character(ranks$value))
        weights <- NULL
        for (i in 1:nrow(ranks))
        {
          if (ranks$count[i] == 1) {
            weights <- c(weights, sum(ranks[1:i, "count"])) ## if a singleton, just add up previous number of values to get rank
          } else {
            start <- sum(ranks[1:i - 1, "count"]) ## get starting point (previous rank)
            ranknum <- 0
            for (j in 1:ranks$count[i])
            {
              ranknum <- ranknum + start + j
            }

            ranknum <- ranknum / ranks$count[i]
            weights <- c(weights, ranknum)
          }
        }

        ranks2 <- cbind(ranks, weights)
        names(ranks2)[3] <- "rank"

        ## populate ex with ranks
        ## need to recast as numeric, may crash otherwise (weird!)
        ex[, nval_idx] <- as.character(ex[, nval_idx])
        ex[, nval_idx] <- as.numeric(ex[, nval_idx])

        ex2 <- merge(ex, ranks2, by.x = as.character(unlist(formula[[2]])), by.y = "value", all.x = TRUE)
        dose_idx <- which(names(ex2) == dose_name)
        ex2 <- ex2[order(ex2[, dose_idx]), ]
        formula_text <- sprintf("rank ~ %s", dose_name)
        rankMeans <- summaryBy(as.formula(formula_text), data = ex2, FUN = mean)

        ## enumerate doseCount
        count <- 0
        for (i in 1:nrow(rankMeans))
        {
          rankMeans$doseCount[i] <- count
          count <- count + 1
        }


        ## find variance
        V <- (nrow(ex) * (nrow(ex) + 1)) / 12

        ## get crit values ... warnings are for CONTROL group, which is fine
        prob05 <- qnorm(1 - (.05 / (2 * max(rankMeans$doseCount))))
        prob01 <- qnorm(1 - (.01 / (2 * max(rankMeans$doseCount))))

        ## loop through doses, determine significance
        rankMeans$DUNSIGN <- 0
        rankMeans$num <- 0

        for (i in 2:nrow(rankMeans)) ## skip CONTROL (first line)
        {

          ## populate DUNSIGN (flag for mean of treatment group being greater than control group)
          if (rankMeans$rank.mean[i] >= rankMeans$rank.mean[1]) {
            rankMeans$DUNSIGN[i] <- 0
          } else {
            rankMeans$DUNSIGN[i] <- -1
          }

          ## find significance level
          rm_idx <- which(dose_name == names(rankMeans))
          ex_idx <- which(dose_name == names(ex2))
          rankdiff <- abs(rankMeans$rank.mean[i] - rankMeans$rank.mean[1])
          num <- nrow(ex2[ex2[, ex_idx] == rankMeans[i, rm_idx], ])
          comp2 <- V * ((1 / num) + 1 / nrow(ex2[ex2[, ex_idx] == rankMeans[1, rm_idx], ]))
          comp2 <- (comp2 * (1 - correction / (nrow(ex2)^3 - nrow(ex2))))^.5

          sig05 <- rankdiff - (prob05 * comp2)

          ptest <- min(1, (1 - pnorm(rankdiff / comp2)) * (2 * max(rankMeans$doseCount)))
          #  message(sprintf("%f %f %f %f %f",rankdiff,comp2,prob05,sig05,ptest))
          sig01 <- rankdiff - (prob01 * comp2)
          signific <- 0
          if (sig05 > 0) {
            signific <- 1
          }

          if (sig01 > 0) {
            signific <- 2
          }

          rankMeans$num[i] <- num
          rankMeans$pvalue[i] <- ptest
        }

        results <- rankMeans

        results_len <- ncol(results)
        temp_names <- colnames(dunn)[unlist(colnames(dunn)) %in% unlist(colnames(data))]

        for (ii in 1:length(temp_names)) {
          idx <- which(temp_names[ii] == colnames(ex2))
          results[, ii + results_len] <- ex2[1, idx]
        }
        A <- results[, results_len + 1:length(temp_names)]
        colnames(A) <- temp_names
        results <- cbind(results[, 1:results_len], A)

        dunn_results <- rbind(dunn_results, results)
      }
    }
  }

  ## check for existence
  if (!is.null(dunn_results)) {
    ## coerce into same format as SHIRLEY results
    dunn_results$TEST <- "DUNN"

    idx <- which(colnames(dunn_results) == dose_name)
    ## cut out controls
    dunn_results <- dunn_results[as.numeric(dunn_results[, idx]) != 0, ]


    names_to_drop <- colnames(dunn_results)
    temp <- sprintf("%sCount", dose_name)
    dose_idx <- which(colnames(dunn_results) == dose_name)
    p_value_idx <- which(colnames(dunn_results) == "pvalue")
    test_idx <- which(colnames(dunn_results) == "TEST")
    remain_idx <- which(!(1:ncol(dunn_results) %in% c(dose_idx, p_value_idx, test_idx)))
    dunn_results <- dunn_results[, c(dose_idx, remain_idx, test_idx, p_value_idx)]
    dunn_results <- dunn_results[, -which(names_to_drop %in% c("rank.mean", temp, "DUNSIGN", "num"))]
  }
  class(dunn_results) <- "ntp.dunn"
  return(dunn_results)
}

## ------------------------
## DUNNETT'S TEST
## ------------------------
#' @title ntp_dunett Dunnett's test
#' @param formula An equation of the form \eqn{Y \sim X.} Here the variable
#' \eqn{Y} is the response of interest, and \eqn{X} represents discrete experimental
#' conditions. For example, if weight is the dependent variable, and you are
#' interested in looking at the trend across sex one would have 'weight ~ sex'.
#' @param data A data frame with column names in the formula.
#' @param dose_name The name of the variable containing the doses in the data frame \eqn{data}.
#' It is expected multiple doses for each of the experimental conditions \eqn{X}.
#' @return The results of Dunnet's  test for each level in \eqn{dose_name}
#' @examples
#' a <- ntp_dunnett(
#'        response ~ sex + response_type,
#'        data = ntp_599_hemotology, 
#'        dose_name = "concentration"
#' )
#' summary(a)
ntp_dunnett <- function(formula, data, dose_name = "dose") {
  mult_comp_test <- NULL
  data[, c(dose_name)] <- as.numeric(unlist(data[, c(dose_name)]))
  temp_str <- strsplit(as.character(formula)[3], " ")[[1]]
  temp_str <- temp_str[temp_str != "+"]
  data <- data[order(data[, c(dose_name)]), ]

  for (ii in 1:length(temp_str)) {
    data <- data[order(data[, temp_str[ii]]), ]
  }

  if (!(dose_name %in% colnames(data))) {
    stop(sprintf("Dose name %s does not appear in the data frame.", dose_name))
  }

  jonck_data <- ntp_jonckeere(formula,
    dose_name = dose_name,
    data = data, pair = "Williams"
  )


  ## loop through all groups flagged as DUNNETT in jonck
  dunnett <- subset(jonck_data, mult_comp_test == "DUNNETT")
  dunnett_results <- NULL


  temp_colnames <- unlist(c(unlist(colnames(dunnett)), dose_name, as.character(unlist(formula[[2]]))))
  temp <- colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[, temp == TRUE])

  if (nrow(dunnett) == 0) {
    warning("The Jonckherre test did not suggest the test be performed. Returning a NULL dataset.")
  }

  if (nrow(dunnett) > 0) {
    for (d in 1:nrow(dunnett)) {
      temp_names <- rep(NA, ncol(dunnett))
      for (j in 1:length(temp_names)) {
        temp_names[j] <- unlist(dunnett[d, j])
      }
      ## KRS - changed "phase" to "phasetype"
      temp_dd <- temp_d
      temp_dd[is.na(temp_dd)] <- ""
      # subset the data based upon the columns that exist in the formula
      for (j in 1:length(dunnett)) {
        myt <- which(colnames(temp_d) == colnames(dunnett)[j])

        if (length(myt) > 0) {
          temp_dd <- temp_dd[as.character(temp_dd[, myt]) == as.character(temp_names[j]), ]
        }
      }

      datatemp <- temp_dd
      temp_idx <- which(colnames(datatemp) == dose_name)
      datatemp <- datatemp[order(datatemp[, temp_idx]), ]

      ## duplicate dose classes from old SAS conversion
      datatemp$dose3 <- as.factor(datatemp[, temp_idx])

      if (length(unique(datatemp[, temp_idx])) > 1) {
        temp_idx2 <- which(!(colnames(dunnett) %in% c("tau", "pvalue", "mult_comp_test")))
        trts <- unique(datatemp[, temp_idx])[-1] ## get unique treatment doses, this will be used in inner loop below
        formula_text <- sprintf("%s ~ dose3", as.character(formula[[2]]))
        data.aov <- aov(as.formula(formula_text), data = datatemp)
        data.dunn <- glht(data.aov, linfct = mcp(dose3 = "Dunnett"))

        ## loop through doses within sex/endpoint combo
        for (j in 1:length(summary(data.dunn)$test$pvalues)) {
          line <- c(temp_names[temp_idx2], trts[j], summary(data.dunn)$test$tstat[j], summary(data.dunn)$test$pvalues[j])
          dunnett_results <- rbind(dunnett_results, line)
        }
        colnames(dunnett_results) <- c(colnames(dunnett)[temp_idx2], dose_name, "tstat", "pvalue")
      }
    }

    if (!is.null(dunnett_results)) {
      rownames(dunnett_results) <- NULL
      dunnett_results <- as.data.frame(dunnett_results)

      dunnett_results$pvalue <- as.numeric(as.character(dunnett_results$pvalue))
      dunnett_results$mult_comp_test <- "DUNNETT"

      ## determine how many asterisks each row deserves
      dunnett_results$mult_comp_signif <- ifelse(dunnett_results$pvalue <= .01, 2,
        ifelse(dunnett_results$pvalue <= .05, 1, 0)
      )


      temp_idx3 <- which(colnames(dunnett_results) == dose_name)
      dunnett_results[, temp_idx3] <- as.numeric(dunnett_results[, temp_idx3])
    }
  }
  class(dunnett_results) <- "ntp.dunnett"
  return(dunnett_results)
}

## ----------------------
## SHIRLEY'S TEST
## ----------------------
#' @title ntp_shirley Shirley's test as programmed at the NTP
#' @param formula An equation of the form \eqn{Y \sim X.} Here the variable
#' \eqn{Y} is the response of interest, and \eqn{X} represents discrete experimental
#' conditions. For example, if weight is the dependent variable, and you are
#' interested in looking at the trend across sex one would have 'weight ~ sex'.
#' @param data A data frame with column names in the formula.
#' @param dose_name The name of the variable containing the doses in the data frame \eqn{data}.
#' It is expected multiple doses for each of the experimental conditions \eqn{X}.
#' @return The results of a non-parametric Shirley's isotone test for trend on
#' each level in \eqn{dose_name}. For more information see: \doi{10.2307/2529789}
#' The returned list contains:
#' \itemize{
#'  \item \code{X}: this represents all the class objects on the right hand side of \eqn{ Y \sim X} above.
#'  \item \code{dose:} the dose groups relative to control.
#'  \item \code{testStats}: Value of the Shirley test statistic.
#'  \item \code{mult_comp_signif}: Test's significance as 0, 1, or 2 which is not-significant,
#'  significant at the 0.05% level and significant at the 0.01% level.
#'  \item \code{mult_comp_test}: The type of test, i.e. "SHIRLEY"
#' }
#' @examples
#' a <- ntp_shirley(weight ~ sex, data = ntp_weight_data)
#' summary(a)
ntp_shirley <- function(formula, data, dose_name = "dose") {
  mult_comp_test <- numTies <- NULL
  data[, c(dose_name)] <- as.numeric(unlist(data[, c(dose_name)]))
  temp_str <- strsplit(as.character(formula)[3], " ")[[1]]
  temp_str <- temp_str[temp_str != "+"]
  data <- data[order(data[, c(dose_name)]), ]

  for (ii in 1:length(temp_str)) {
    data <- data[order(data[, temp_str[ii]]), ]
  }

  if (!(dose_name %in% colnames(data))) {
    stop(sprintf("Dose name %s does not appear in the data frame.", dose_name))
  }

  ## loop through all groups flagged as SHIRLEY in jonck
  jonck_data <- ntp_jonckeere(formula,
    dose_name = dose_name,
    data = data, pair = "Shirley"
  )
  shirley_results <- NULL
  shirley <- subset(jonck_data, mult_comp_test == "SHIRLEY")

  temp_colnames <- unlist(c(unlist(colnames(shirley)), dose_name, as.character(unlist(formula[[2]]))))
  temp <- colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[, temp == TRUE])

  if (nrow(shirley) > 0) {
    for (s in 1:nrow(shirley))
    {
      testStats <- NULL
      doseCount <- NULL
      dose <- NULL
      num <- NULL

      temp_names <- rep(NA, ncol(shirley))
      for (j in 1:length(temp_names)) {
        temp_names[j] <- unlist(shirley[s, j])
      }


      ## KRS - changed "phase" to "phasetype"
      temp_dd <- temp_d
      temp_dd[is.na(temp_dd)] <- ""
      # subset the data based upon the columns that exist in the formula
      for (j in 1:length(shirley)) {
        myt <- which(colnames(temp_d) == colnames(shirley)[j])

        if (length(myt) > 0) {
          temp_dd <- temp_dd[as.character(temp_dd[, myt]) == as.character(temp_names[j]), ]
        }
      }
      ex <- temp_dd
      dose_idx <- which(colnames(ex) == dose_name)
      nval_idx <- which(colnames(ex) == as.character(unlist(formula[[2]])))
      ex[, dose_idx] <- as.numeric(ex[, dose_idx])

      ## make sure there are multiple doses to compare, and control is present
      if ((length(unique(ex[, dose_idx])) > 1) & (0 %in% unique(ex[, dose_idx]))) {
        exLoop <- ex

        for (g in (length(unique(ex[, dose_idx])) - 1):1) {
          ## get ties for use in correction formula
          ties <- as.data.frame(table(table(exLoop[, nval_idx])))

          if (nrow(ties) > 1) {
            names(ties) <- c("numTies", "count")
            ties$numTies <- as.integer(as.character(ties$numTies))

            ## remove singletons
            ties <- subset(ties, numTies != 1)

            if (nrow(ties) > 0) {
              correction <- 0
              for (i in 1:nrow(ties))
              {
                correction <- correction + (ties$count[i] * (ties$numTies[i]^3 - ties$numTies[i]))
              }

              correction <- correction / (12 * (nrow(exLoop) - 1))
            } else {
              correction <- 0
            }
          } else {
            correction <- 0
          }

          ## rankings ... lowest value gets rank of 1
          ranks <- as.data.frame(table(exLoop[, nval_idx]))
          names(ranks) <- c("value", "count")
          ranks$value <- as.numeric(as.character(ranks$value))

          weights <- NULL
          for (i in 1:nrow(ranks)) {
            if (ranks$count[i] == 1) {
              weights <- c(weights, sum(ranks[1:i, "count"])) ## if a singleton, just add up previous number of values to get rank
            } else {
              start <- sum(ranks[1:i - 1, "count"]) ## get starting point (previous rank)
              ranknum <- 0
              for (j in 1:ranks$count[i])
              {
                ranknum <- ranknum + start + j
              }

              ranknum <- ranknum / ranks$count[i]
              weights <- c(weights, ranknum)
            }
          }

          ranks2 <- cbind(ranks, weights)
          names(ranks2) <- c("value", "count", "rank")

          ## populate exLoop with ranks
          ## need to recast as numeric, may crash otherwise (weird!)
          exLoop[, nval_idx] <- as.character(exLoop[, nval_idx])
          exLoop[, nval_idx] <- as.numeric(exLoop[, nval_idx])

          exLoop2 <- merge(exLoop, ranks2, by.x = as.character(unlist(formula[[2]])), by.y = "value", all.x = TRUE)
          exLoop2 <- exLoop2[order(exLoop2[, dose_idx]), ]

          formula_text <- sprintf("rank ~ %s", dose_name)

          rankSums <- summaryBy(as.formula(formula_text), data = exLoop2, FUN = c(sum, length)) ## !!
          rankMeans <- summaryBy(as.formula(formula_text), data = exLoop2, FUN = c(mean, length)) ## raw means, will overwrite with corrected means after next step

          mean_line <- NULL
          numer <- 0
          denom <- 0
          for (i in nrow(rankSums):1) {
            numer <- numer + rankSums$rank.sum[i]
            denom <- denom + rankSums$rank.length[i]

            if (i == 1) {
              mean_temp <- rankSums$rank.sum[i] / rankSums$rank.length[i]
            } else {
              mean_temp <- numer / denom
            }

            mean_line <- c(mean_line, mean_temp)
          }

          rankMeans$rank.mean <- rev(mean_line)
          ## find test statistic
          V <- (nrow(exLoop) * (nrow(exLoop) + 1)) / 12 - correction
          Ri <- nrow(exLoop[exLoop[, which(colnames(exLoop) == dose_name)] == rankMeans[nrow(rankMeans), which(colnames(rankMeans) == dose_name)], ])
          C <- nrow(exLoop[exLoop[, which(colnames(exLoop) == dose_name)] == rankMeans[1, which(colnames(exLoop) == dose_name)], ])

          if (shirley$tau[s] >= 0) {
            dosemean <- max(rankMeans$rank.mean[2:nrow(rankMeans)])
            shrl_num <- dosemean - rankMeans$rank.mean[1]
          } else {
            dosemean <- min(rankMeans$rank.mean[2:nrow(rankMeans)])
            shrl_num <- rankMeans$rank.mean[1] - dosemean
          }

          T <- shrl_num * (V * (1 / Ri + 1 / C))^-.5 ## shirlstat in SAS code

          testStats <- c(testStats, T)
          doseCount <- c(doseCount, g)

          dose <- c(dose, unique(ex[, which(colnames(ex) == dose_name)])[g + 1])
          num <- c(num, sum(ex[, which(colnames(ex) == dose_name)] == unique(ex[, which(colnames(ex) == dose_name)])[g + 1]))

          td <- which(colnames(exLoop) == dose_name)
          exLoop <- exLoop[exLoop[, td] != unique(exLoop[, td])[length(unique(exLoop[, td]))], ]
          ## remove latest dosage before returning to top of loop
        }


        tshirl <- shirley[, -which(colnames(shirley) %in% c("tau", "pvalue", "mult_comp_test")), drop = F][s, , drop = F]

        results <- as.data.frame(cbind(dose, num, doseCount, testStats))

        real_temp <- as.data.frame(matrix(NA, nrow = nrow(results), ncol = ncol(tshirl)))

        for (ii in 1:ncol(tshirl)) {
          real_temp[, ii] <- tshirl[ii]
        }

        names(real_temp) <- colnames(tshirl)
        results <- cbind(real_temp, results)


        ## add SAS crit values
        C01 <- c(0, 2.575, 2.607, 2.615, 2.618, 2.620, 2.621, 2.622)
        C05 <- c(0, 1.96, 2.015, 2.032, 2.040, 2.044, 2.047, 2.0485)

        B01 <- c(0, 0, 3, 4, 4, 4, 4, 4)
        B05 <- c(0, 0, 3, 4, 5, 6, 6, 6)


        ## find crit values, generate number of stars
        results$mult_comp_signif <- 0
        nonsignif_flag <- "NO" ## to match Laura's version force all doses after first non-signif dose to zero

        if (length(doseCount) <= 7) { ## cannot handle more than 7 non control groups ... no crit values to compare to
          for (i in 1:nrow(results))
          {
            if (nonsignif_flag == "NO") {
              dosenum <- results$doseCount[i] + 1
              crit05 <- C05[dosenum] - ((B05[dosenum] / 100) * (1 - (results$num[i] / results$num[nrow(results)])))
              crit01 <- C01[dosenum] - ((B01[dosenum] / 100) * (1 - (results$num[i] / results$num[nrow(results)])))
              if (results$testStats[i] >= crit01) {
                results$mult_comp_signif[i] <- 2
              } else {
                if (results$testStats[i] >= crit05) {
                  results$mult_comp_signif[i] <- 1
                } else {
                  results$mult_comp_signif[i] <- 0
                  nonsignif_flag <- "YES"
                }
              }
            }
          }
        } else {
          # results$mult_comp_signif <- 'NA'
          results$mult_comp_signif <- ""
        }

        if (!is.null(shirley_results)) {
          shirley_results <- rbind(shirley_results, results)
        } else {
          shirley_results <- results
        }
      }
    }
  }
  ## check for existence

  if (!is.null(shirley_results)) {


    ## coerce into same format as SHIRLEY results
    shirley_results$mult_comp_test <- "SHIRLEY"
    idx <- c(which(colnames(shirley_results) %in% colnames(tshirl)), which(colnames(shirley_results) == "dose"))
    # way too complicated loop to do something simple *sigh*
    for (ii in length(idx):1) {
      reorder <- sort(shirley_results[, idx[ii]], index = TRUE)$ix
      shirley_results <- shirley_results[reorder, ]
    }
    shirley_results <- shirley_results[, -which(colnames(shirley_results) %in% c("num", "doseCount"))]
    ## remove NA's
    shirley_results[is.na(shirley_results)] <- ""
  }

  class(shirley_results) <- "ntp.shirley"

  return(shirley_results)
}

.summary_ntpwilliams <- function(object, ...) {
  class(object) <- "data.frame"
  cat("Williams Trend Test: Monotone Change from control?\n")
  cat("--------------------------------------------------\n")
  loc <- which(names(object) %in% c("willStat", "mult_comp_signif", "mult_comp_test"))
  data_one <- object[, -loc]
  data_two <- object[, c("mult_comp_signif")]
  data_a <- rep("No", length(data_two))
  data_a[data_two == 1] <- "<0.05"
  data_a[data_two == 2] <- "<0.01"
  output <- data.frame(data_one, data_a)
  names(output) <- c(names(data_one), "Significant")

  print(output, row.names = F)
}

.summary_ntpdunn <- function(object, ...) {
  class(object) <- "data.frame"
  loc <- which(names(object) == "TEST")
  object <- object[, -loc]
  pv_loc <- which(names(object) == "pvalue")
  data_two <- object[, pv_loc]
  data_one <- object[, -c(loc, pv_loc)]


  cat("Dunn Trend Test: Significant Change from control? \n")
  cat("--------------------------------------------------\n")

  data_a <- rep("No", length(data_two))
  data_a[data_two < 0.05] <- "<0.05"
  data_a[data_two < 0.01] <- "<0.01"
  output <- data.frame(data_one, data_a)
  names(output) <- c(names(data_one), "Significant")
  print(output, row.names = FALSE)
}

.summary_ntpdunnett <- function(object, ...) {
  class(object) <- "data.frame"
  loc <- which(names(object) %in% c("TEST", "tstat", "mult_comp_test", "tstat", "mult_comp_signif"))
  object <- object[, -loc]
  pv_loc <- which(names(object) == "pvalue")
  data_two <- object[, pv_loc]
  data_one <- object[, -c(loc, pv_loc)]

  cat("Dunnett Trend Test: Significant Change from control? \n")
  cat("---------------------------------------------------  \n")

  data_a <- rep("No", length(data_two))
  data_a[data_two < 0.10] <- "<0.10"
  data_a[data_two < 0.05] <- "<0.05"
  data_a[data_two < 0.01] <- "<0.01"
  output <- data.frame(data_one, data_a)
  names(output) <- c(names(data_one), "Significant")
  print(output, row.names = FALSE)
}

.summary_ntpshirley <- function(object, ...) {
  class(object) <- "data.frame"
  loc <- which(names(object) %in% c("testStats", "mult_comp_test"))
  object <- object[, -loc]
  pv_loc <- which(names(object) == "mult_comp_signif")
  data_two <- object[, pv_loc]
  data_one <- object[, -c(loc, pv_loc)]

  cat("Shriley's Trend Test: Monotone Change from control?  \n")
  cat("---------------------------------------------------  \n")

  data_a <- rep("No", length(data_two))
  data_a[data_two == 1] <- "<0.05"
  data_a[data_two == 2] <- "<0.01"
  output <- data.frame(data_one, data_a)
  names(output) <- c(names(data_one), "Significant")
  print(output, row.names = FALSE)
}
