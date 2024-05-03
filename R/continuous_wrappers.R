#/*
# Copyright 2020  NIEHS <matt.wheeler@nih.gov>
#*Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
#*and associated documentation files (the "Software"), to deal in the Software without restriction, 
#*including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#*sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
#*is furnished to do so, subject to the following conditions:
# *
#*The above copyright notice and this permission notice shall be included in all copies 
#*or substantial portions of the Software.

#*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#*INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#*PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#*HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#*CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#*OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#' Fit a single continuous BMD model.
#'
#' @title single_continuous_fit - Fit a single continuous BMD model.
#' @param D doses matrix
#' @param Y response matrix
#' @param model_type Specifies the mean model, or log-mean 
#' for lognormal data.  For more information, see Aerts, Wheeler, and Abrahantes (2021)
#'
#'  \itemize{
#'      \item \code{"exp-aerts"}:         \eqn{f(x) = a(1 + (c-1)(1-\exp(-bx^{d}))) }
#'      \item \code{"invexp-aerts"}:      \eqn{f(x) = a(1 + (c-1)(\exp(-bx^{-d})))}
#'      \item \code{"hill-aerts"}:        \eqn{f(x) = a(1 + (c-1)(1-\frac{b^d}{b^d + x^d}))}
#'      \item \code{"gamma-aerts"}:       \eqn{f(x) = a(1 + (c-1)(\Gamma(bx^d;\xi)))}
#'      \item \code{"invgamma-aerts"}:    \eqn{f(x) = a(1 + (c-1)(1-\Gamma(bx^{-d};\xi)))}
#'      \item \code{"lomax-aerts"}:       \eqn{f(x) = a\left\{1 + (c-1)(1-\left(\frac{b}{b+x^d} \right)^\xi) \right\}}
#'      \item \code{"invlomax-aerts"}:    \eqn{f(x) = a\left\{1 + (c-1)(\left(\frac{b}{b+x^{-d}} \right))^\xi \right\}}
#'      \item \code{"lognormal-aerts"}:   \eqn{f(x) = a\left\{1 + (c-1)\left(\Phi( \ln(b) + d\times \ln(x))\right) \right\}}
#'      \item \code{"logskew-aerts"}:     \eqn{f(x) = a\left\{1 + (c-1)\left(\Phi_{SN}( \ln(b) + d\times \ln(x); \xi )\right) \right\}}
#'      \item \code{"invlogskew-aerts"}:  \eqn{f(x) = a\left\{1 + (c-1)\left(1 - \Phi_{SN}( \ln(b) - d\times \ln(x); \xi )\right) \right\}}
#'      \item \code{"logistic-aerts"}:    \eqn{f(x) = \frac{c}{1 + \exp(-a - b\times x^d)} }
#'      \item \code{"probit-aerts"}:      \eqn{f(x) = c\left(\Phi(a + b\times x^d)\right) }
#'      \item \code{"gamma-efsa"}:        \eqn{f(x) = a(1 + (c-1)(\Gamma(bx; d))) }
#'      \item \code{"LMS"}:               \eqn{f(x) = a(1 + (c-1)(1 - \exp(-bx - dx^2))) }
#'    }
#'   Here: \eqn{\Phi(\cdot)} is the standard normal distribution and
#'         \eqn{\Phi_{SN}(\cdot;\cdot)} is the skew-normal distribution
#'
#'    Legacy continuous models based upon US EPA BMDS software. 
#'    These models can be fit using Bayesian and Maximum Likelihood. 
#'    \itemize{
#'      \item \code{"hill"}: \eqn{f(x) = a + b\frac{x^d}{c^d + x^d}}
#'      \item \code{"exp-3"}: \eqn{f(x) = a(1-exp(-bx^{d})) }
#'      \item \code{"exp-5"}: \eqn{f(x) = a(1 + (c-1)(1-exp(-bx^{d}))) }
#'      \item \code{"power"}: \eqn{f(x) = a + b x^d}
#'      \item \code{"polynomial"}: \eqn{\beta_0 + \beta_1 x + \ldots \beta_p x^p}, where p is equal to 
#'      \code{degree} (defined) below. 
#'    }
#' @param fit_type the method used to fit (laplace, mle, or mcmc)
#' @param prior Prior / model for the continuous fit. If this is specified, it overrides the parameters 'model_type' and 'distribution.' 
#' @param BMR_TYPE Specifies the type of benchmark dose analysis to be performed. For continuous models, there are four types of BMD definitions that are commonly used:
#'    \itemize{
#'      \item Standard deviation is the default option, but it can be explicitly specified with 'BMR_TYPE = "sd"' This definition defines the BMD as the dose associated with the mean/median changing a specified number of standard deviations from the mean at the control dose., i.e., it is the dose, BMD, that solves \eqn{\mid f(dose)-f(0) \mid = BMR \times \sigma}
#'      \item Relative deviation can be specified with 'BMR_TYPE = "rel"'. This defines the BMD as the dose that changes the control mean/median a certain percentage from the background dose, i.e. it is the dose, BMD that solves \eqn{\mid f(dose) - f(0) \mid = (1 \pm BMR) f(0)} 
#'      \item Hybrid deviation can be specified with 'BMR_TYPE = "hybrid"'.  This defines the BMD that changes the probability of an adverse event by a stated amount relitive to no exposure (i.e 0).  That is, it is the dose, BMD, that solves \eqn{\frac{Pr(X > x| dose) - Pr(X >x|0)}{Pr(X < x|0)} = BMR}. For this definition, \eqn{Pr(X < x|0) = 1 - Pr(X > X|0) = \pi_0}, where \eqn{0 \leq \pi_0 < 1} is defined by the user as "point_p," and it defaults to 0.01.  Note: this discussion assumed increasing data.  The fitter determines the direction of the data and inverts the probability statements for decreasing data.
#'      \item Absolute deviation can be specified with 'BMR_TYPE="abs"'. This defines the BMD as an absolute change from the control dose of zero by a specified amount. That is the BMD is the dose that solves the equation \eqn{\mid f(dose) - f(0) \mid = BMR}. 
#'    }  
#' @param BMR This option specifies the benchmark response BMR. The BMR is defined in relation to the BMD calculation requested (see BMD).  By default, the "BMR = 0.1."\cr
#' @param point_p This option is only used for hybrid BMD calculations. It defines a probability that is the cutpoint for observations.  It is the probability that observations have this probability, or less, of being observed at the background dose. \cr
#' @param alpha Alpha is the specified nominal coverage rate for computation of the lower bound on the BMDL and BMDU, i.e., one computes a \eqn{100\times(1-\alpha)\%} confidence interval.  For the interval (BMDL,BMDU) this is a \eqn{100\times(1-2\alpha)\%} confidence interval.  By default, it is set to 0.05.
#' @param samples the number of samples to take (MCMC only)
#' @param degree the number of degrees of a polynomial model. Only used for polynomial models. 
#' @param burnin the number of burnin samples to take (MCMC only)
#' @param distribution The underlying distribution used as the data distribution. 
#' @param BMD_priors Logical specifying if BMD Express-3 priors should be used instead of ToxicR defaults
#' @param ewald perform Wald CI computation instead of the default profile likelihood computation. This is the the 'FAST BMD' method of Ewald et al (2021)
#' @param transform Transforms doses using \eqn{\log(dose+\sqrt{dose^2+1})}. Note: this is a log transform that has a derivative defined when dose =0.
#' @param BMD_TYPE Deprecated version of BMR_TYPE that specifies the type of benchmark dose analysis to be performed
#' @param threads specify the number of OpenMP threads to use for the calculations. Default = 2
#' @param seed specify the GSL seed. Default = 12331
#' @return Returns a model object class with the following structure:
#' \itemize{
#'    \item \code{full_model}:  The model along with the likelihood distribution. 
#'    \item \code{bmd}:  A vector containing the benchmark dose (BMD) and \eqn{100\times(1-2\alpha)} confidence intervals. 
#'    \item \code{parameters}: The parameter estimates produced by the procedure, which are relative to the model '
#'                             given in \code{full_model}.  The last parameter is always the estimate for \eqn{\log(sigma^2)}.
#'    \item \code{covariance}: The variance-covariance matrix for the parameters.  
#'    \item \code{bmd_dis}:  Quantiles for the BMD distribution. 
#'    \item \code{maximum}:  The maximum value of the likelihod/posterior. 
#'    \item \code{Deviance}:  An array used to compute the analysis of deviance table. 
#'    \item \code{prior}:     This value gives the prior for the Bayesian analysis. 
#'    \item \code{model}:     Parameter specifies t mean model used. 
#'    \item \code{options}:   Options used in the fitting procedure.
#'    \item \code{data}:     The data used in the fit. 
#'    \item \code{transformed}: Are the data \eqn{\log(x+\sqrt{x^2+1})} transformed? 
#'    \itemize{
#'        When MCMC is specified, an additional variable \code{mcmc_result} 
#'        has the following two variables:
#'        \item \code{PARM_samples}:  matrix of parameter samples. 
#'        \item \code{BMD_samples}: vector of BMD sampled values. 
#'    }
#' }
#' 
#' @examples 
#' M2           <- matrix(0,nrow=5,ncol=4)
#' colnames(M2) <- c("Dose","Resp","N","StDev")
#' M2[,1] <- c(0,25,50,100,200)
#' M2[,2] <- c(6,5.2,2.4,1.1,0.75)
#' M2[,3] <- c(20,20,19,20,20)
#' M2[,4] <- c(1.2,1.1,0.81,0.74,0.66)
#' model = single_continuous_fit(
#'          M2[,1,drop=FALSE], 
#'          M2[,2:4], 
#'          BMR_TYPE="sd", 
#'          BMR=1, ewald = TRUE,
#'          distribution = "normal",
#'          fit_type="laplace",
#'          model_type = "hill",
#'          threads = 2, 
#'          seed = 12331
#' )
#' 
#' summary(model)
single_continuous_fit <- function(D,Y,model_type="hill", fit_type = "laplace",
                                   prior=NA, BMR_TYPE = "sd", 
                                   BMR = 0.1, point_p = 0.01, distribution = "normal-ncv",
                                   alpha = 0.05, samples = 25000, degree=2,
                                   burnin = 1000, BMD_priors = FALSE, ewald = FALSE,
                                   transform = FALSE, BMD_TYPE = NA, threads = 2, seed = 12331){
    .setseedGSL(seed)
    Y <- as.matrix(Y) 
    D <- as.matrix(D) 
    
    dis_type = which(distribution  == c("normal","normal-ncv","lognormal"))
    
    if (dis_type == 3){
       is_neg = .check_negative_response(Y)
       if (is_neg){
         stop("Can't fit a negative response to the log-normal distribution.")
       }
    }
    
    DATA <- cbind(D,Y);
    test <-  .check_for_na(DATA)
    Y = Y[test==TRUE,,drop=F]
    D = D[test==TRUE,,drop=F]
    DATA <- cbind(D,Y);
   
    
    myD = Y; 
    sstat = F # set sufficient statistics to false if there is only one column
    if (ncol(Y) > 1){
        sstat=T
    }
    
    if (!is.na(prior[1])){
      #parse the prior
      if (!( "BMD_Bayes_continuous_model" %in% class(prior))){
        stop("Prior is not the correct form. Please use a Bayesian Continuous Prior Model.")
      }
      t_prior_result <- .parse_prior(prior)
      distribution <- t_prior_result$distribution
      model_type   <- t_prior_result$model
      prior = t_prior_result$prior
      PR = t_prior_result$prior
    }else{
      dmodel = which(model_type==.continuous_models)
      
      if (identical(dmodel, integer(0))){
        stop('Please specify one of the following model types: \n
            "hill","exp-3","exp-5","power","polynomial", "exp-aerts", "invexp-aerts", "gamma-aerts", "invgamma-aerts", "hill-aerts",
            "lomax-aerts", "invlomax-aerts", "lognormal-aerts", "logskew-aerts", "invlogskew-aerts", "logistic-aerts", "probit-aerts", "LMS",
            "gamma-efsa"')
      }

      PR    = .bayesian_prior_continuous_default(model_type,distribution,degree)
      #specify variance of last parameter to variance of response
      if(distribution == "lognormal"){
           if (ncol(Y)>1){
               PR$priors[nrow(PR$priors),2]= log(mean(Y[,3])^2)
           }else{
               PR$priors[nrow(PR$priors),2] = log(var(log(Y[,1])))
           }
      }else{
           if (ncol(Y)>1){
                if (distribution == "normal"){
                     PR$priors[nrow(PR$priors),2]   = log(mean(Y[,3])^2)
                }else{
                     PR$priors[nrow(PR$priors),2]   = log(abs(mean(Y[1,]))/mean(Y[,3])^2)
                }
           }else{ 
                  if (distribution == "normal"){
                    PR$priors[nrow(PR$priors),2]   = log(var(Y[,1]))
                  }else{
                    PR$priors[nrow(PR$priors),2]   = log(abs(mean(Y[,1]))/var(Y[,1]))
                  }
           }
      }
      t_prior_result = create_continuous_prior(PR,model_type,distribution,degree)
      PR = t_prior_result$prior
    }
    dmodel = which(model_type==.continuous_models)

    fit_type = tolower(fit_type)
    type_of_fit = which(fit_type == c('laplace','mle','mcmc'))
    if (identical(type_of_fit,integer(0))){
      stop("Please choose one of the following fit types: 'laplace','mle','mcmc.' ")
    }

    if(!is.na(BMD_TYPE)){
      warning("BMD_TYPE is deprecated. Please use BMR_TYPE instead")
    }else{
      BMD_TYPE = BMR_TYPE
    }
    
    rt = which(BMD_TYPE==c('abs','sd','rel','hybrid'))
    
    if (rt == 4){
      rt = 6; 
    }

    
    if(identical(dis_type,integer(0))){
      stop('Please specify the distribution as one of the following:\n
            "normal","normal-ncv","lognormal"')
    }


    if (ncol(DATA)==4){
      colnames(DATA) =  c("Dose","Resp","N","StDev")
    }else if (ncol(DATA) == 2){
      colnames(DATA) =  c("Dose","Resp")
    }else{
      stop("The data do not appear to be in the correct format.")
    }
    #permute the matrix to the internal C values
    # Hill = 6, Exp3 = 3, Exp5 = 5, Power = 8, 
    # FUNL = 10, aerts = 11-22
    permuteMat = cbind(c(1,2,3,4,5,6,7:20),c(6,3,5,8,10,666, 11:24))
    fitmodel = permuteMat[dmodel,2]
    if (fitmodel == 10 && dis_type == 3){
         stop("The FUNL model is currently not defined for Log-Normal distribution.")
    }
    if (fitmodel == 666 && dis_type == 3){
         stop("Polynomial models are currently not defined for the Log-Normal distribution.")
    }
    if (fitmodel %in% 11:24 && fit_type == "mle"){
      stop("Aerts models are currently not supported with the frequentist approach.")
    }
    if(any(PR[,1] >= 4)){
      warning("Gamma and PERT priors are still under development and may be incorrect.")
    }
    
    #Fit to determine direction. 
    #This is done for the BMD estimate. 
    model_data = list(); 
    model_data$X = D; model_data$SSTAT = Y
    if (sstat == T){
      temp.fit <- lm(model_data$SSTAT[,1] ~ model_data$X,
                     weights=(1/model_data$SSTAT[,3]^2)*model_data$SSTAT[,2])
    }else{
      temp.fit <- lm(model_data$SSTAT[,1]~model_data$X)
    }
    
    is_increasing = F
    if (coefficients(temp.fit)[2] > 0){
      is_increasing = T
    }
    
    if (!is_increasing){
         if (BMD_TYPE == "rel"){
              BMR = 1-BMR
         }
    }

    #For MLE 
    if (type_of_fit == 2){
      .set_threads(threads)
      PR = .MLE_bounds_continuous(model_type,distribution,degree, is_increasing)
      PR = PR$priors
    }

    if (distribution == "lognormal"){
      is_log_normal = TRUE
    }else{
      is_log_normal = FALSE
    }
    
    if (distribution == "normal-ncv"){
      constVar = FALSE
      vt = 2;
    }else{
      vt = 1
      constVar = TRUE
    }

    if(BMD_priors == TRUE){
      #MLE differences (only for bounds)
      if(type_of_fit == 2){
        if(fitmodel == 8){ #power model has different lower bound on 3rd parameter
          PR[3,4] <- 1
        } 
        if(fitmodel == 666){#polynomial model has many changes
          if(distribution == "normal-ncv"){
            temp <- nrow(PR) - 2 #NCV term
            PR[temp, 5] <- 100 #UB on NCV parameter in 100
            PR[2:(temp-1), 4] <- -10000
            PR[2:(temp-1), 5] <- 10000
          }
          if(distribution == 'normal'){
            PR[1,4:5] <- c(-1000, 1000)
            if(is_increasing){
              PR[2:(nrow(PR)-1), 4] <- 0.00001
              PR[2:(nrow(PR)-1), 5] <- 18
            } else{
              PR[2:(nrow(PR)-1), 5] <- -0.00001
              PR[2:(nrow(PR)-1), 4] <- -18
            }
          }
        }
      } else{
        if(fitmodel == 3){ #exp-3
          PR[2,] <- c(2,0,6.9,0,10000)
        }
        if(fitmodel == 6){ #hill
          PR[2,] <- c(1,0,1000,-10000, 10000)
        }
        if(fitmodel == 5){ #exp-5
          PR[2,] <- c(1,0,1000,-10000, 10000)
          if(distribution == "normal-ncv"){
            PR[5, 3] <- 0.75
            PR[6,4:5] <- c(-18, 18)
          }
        }
        if(fitmodel == 8 & distribution == "normal-ncv"){ #power NCV
          PR[2,] <- c(1,0,10, -10000, 10000)
        }
      }
    }
    

    if (identical(rt,integer(0))){
    	 stop("Please specify one of the following BMRF types:
    		  'abs','sd','rel','hybrid'")
    }
    
    if (rt == 4) {rt = 6;} #internally in Rcpp hybrid is coded as 6	
    
    tmodel <- permuteMat[dmodel,2] 
    
    options <- c(rt,BMR,point_p,alpha, is_increasing,constVar,burnin,samples,transform)
    
    ## pick a distribution type 
    if (is_log_normal){
      dist_type = 3 #log normal
    }else{
      if (constVar){
          dist_type = 1 # normal
      }else{
          dist_type = 2 # normal-ncv
      }
    }
  ##  print(PR)
  # // return(PR)
    if (fit_type == "mcmc"){
      
      .set_threads(threads)
      rvals <- .run_continuous_single_mcmc(fitmodel,model_data$SSTAT,model_data$X,
                                          PR ,options, is_log_normal, sstat) 
   
      if (model_type == "exp-3"){
        rvals$PARMS = rvals$PARMS[,-3]
        rvals$mcmc_result$PARM_samples = rvals$mcmc_result$PARM_samples[,-3]
      }
   
      rvals$bmd      <- c(rvals$fitted_model$bmd,NA,NA) 
      
      rvals$full_model <- rvals$fitted_model$full_model
      rvals$parameters <- rvals$fitted_model$parameters
      rvals$covariance <- rvals$fitted_model$covariance
      rvals$maximum <- rvals$fitted_model$maximum
      
      rvals$prior    <- t_prior_result
      rvals$bmd_dist <- rvals$fitted_model$bmd_dist
      if (!identical(rvals$bmd_dist, numeric(0))){
        temp_me = rvals$bmd_dist
        #rarely gives error, so wrap in tryCatch: if fails, create 2x2 0 matrix
        tryCatch({
          temp_me = temp_me[!is.infinite(temp_me[,1]),]
          temp_me = temp_me[!is.na(temp_me[,1]),]
          temp_me = temp_me[!is.nan(temp_me[,1]),]
        }, error = function(e) temp_me <- matrix(0, nrow = 2, ncol = 2))
        #nrow can throw error too
        if(is.null(nrow(temp_me))){
          temp_me <- matrix(0, nrow = 2, ncol = 2)
        }
        if( nrow(temp_me) > 5){
          te <- splinefun(temp_me[,2],temp_me[,1],method="monoH.FC",ties=mean)
          rvals$bmd[2:3]  <- c(te(alpha),te(1-alpha))
        }else{
          rvals$bmd[2:3] <- c(NA,NA)
        }
      }
      #if less than 10 unique jumps, warn that mixing was poor
      if(length(unique(rvals$mcmc_result$BMD_samples)) < 10){
        warning("MCMC did not mix well!", call. = FALSE)
      }
      
      class(rvals) <- "BMDcont_fit_MCMC"; 
      rvals$model  <- model_type
      rvals$options <- options
      rvals$data <- DATA
      names(rvals$bmd) <- c("BMD","BMDL","BMDU")
  
      rvals$prior <- t_prior_result
      rvals$transformed <- transform
      class(rvals) <- "BMDcont_fit_MCMC"
      rvals$fitted_model <- NULL
      return(rvals)
    }else{
      
      options[7] <- (ewald == TRUE)*1
      .set_threads(threads)
      rvals   <- .run_continuous_single(fitmodel,model_data$SSTAT,model_data$X,
  						                          PR,options, dist_type)
     
      rvals$bmd_dist = rvals$bmd_dist[!is.infinite(rvals$bmd_dist[,1]),,drop=F]
      rvals$bmd_dist = rvals$bmd_dist[!is.na(rvals$bmd_dist[,1]),,drop=F]
      rvals$bmd     <- c(rvals$bmd,NA,NA)
      names(rvals$bmd) <- c("BMD","BMDL","BMDU")
      if (fit_type == "laplace"){
        rvals$prior    <- t_prior_result
      }
      if (!identical(rvals$bmd_dist, numeric(0))){
        temp_me = rvals$bmd_dist
        tryCatch({
          temp_me = temp_me[!is.infinite(temp_me[,1]),]
          temp_me = temp_me[!is.na(temp_me[,1]),]
          temp_me = temp_me[!is.nan(temp_me[,1]),]
        }, error = function(e) temp_me <- matrix(0, nrow = 2, ncol = 2))
        #nrow can throw error too
        if(is.null(nrow(temp_me))){
          temp_me <- matrix(0, nrow = 2, ncol = 2)
        }
        if( nrow(temp_me) > 5){
          te <- splinefun(temp_me[,2],temp_me[,1],method="monoH.FC",ties=mean)
          rvals$bmd[2:3]  <- c(te(alpha),te(1-alpha))
        }else{
          rvals$bmd[2:3] <- c(NA,NA)
        }
      }
      names(rvals$bmd) <- c("BMD","BMDL","BMDU")
      rvals$model   <- model_type
      rvals$options <- options
      rvals$data    <- DATA
      class(rvals)  <- "BMDcont_fit_maximized"
      rvals$transformed <- transform
      return (rvals)
    }
}
