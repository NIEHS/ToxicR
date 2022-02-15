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
#* 
#  * 
#  */

#' Fit a single continuous BMD model.
#'
#' @title single_continuous_fit - Fit a single continuous BMD model.
#' @param D doses matrix
#' @param Y response matrix
#' @param fit_type the method used to fit (laplace, mle, or mcmc)
#' @param BMD_TYPE BMD_TYPE specifies the type of benchmark dose analysis to be performed. For continuous models, there are four types of BMD definitions that are commonly used. \cr
#' -	Standard deviation is the default option, but it can be explicitly specified with 'BMR_TYPE = "sd"' This definition defines the BMD as the dose associated with the mean/median changing a specified number of standard deviations from the mean at the control dose., i.e., it is the dose, BMD, that solves \eqn{\mid f(dose)-f(0) \mid = BMR \times \sigma} \cr
#' -	Relative deviation can be specified with 'BMR_TYPE = "rel"'. This defines the BMD as the dose that changes the control mean/median a certain percentage from the background dose, i.e. it is the dose, BMD that solves \eqn{\mid f(dose) - f(0) \mid = (1 \pm BMR) f(0)} \cr
#' -	Hybrid deviation can be specified with 'BMR_TYPE = "hybrid"'.  This defines the BMD that changes the probability of an adverse event by a stated amount relitive to no exposure (i.e 0).  That is, it is the dose, BMD, that solves \eqn{\frac{Pr(X > x| dose) - Pr(X >x|0)}{Pr(X < x|0)} = BMR}. For this definition, \eqn{Pr(X < x|0) = 1 - Pr(X > X|0) = \pi_0}, where \eqn{0 \leq \pi_0 < 1} is defined by the user as "point_p," and it defaults to 0.01.  Note: this discussion assumed increasing data.  The fitter determines the direction of the data and inverts the probability statements for decreasing data. \cr
#' -	Absolute deviation can be specified with 'BMR_TYPE="abs"'. This defines the BMD as an absolute change from the control dose of zero by a specified amount. That is the BMD is the dose that solves the equation \eqn{\mid f(dose) - f(0) \mid = BMR}  
#' @param BRM This option specifies the benchmark response BMR. The BMR is defined in relation to the BMD calculation requested (see BMD).  By default, the "BMR = 0.1."
#' @param point_p This option is only used for hybrid BMD calculations. It defines a probability that is the cutpoint for observations.  It is the probability that observations have this probability, or less, of being observed at the background dose.
#' @param alpha Alpha is the specified nominal coverage rate for computation of the lower bound on the BMDL and BMDU, i.e., one computes a \eqn{100\times(1-\alpha)% confidence interval}.  For the interval (BMDL,BMDU) this is a \eqn{100\times(1-2\alpha)% confidence interval}.  By default, it is set to 0.05.
#' @param samples the number of samples to take (MCMC only)
#' @param degree the number of degrees of a polynomial model. Only used for polynomial models. 
#' @param burnin the number of burnin samples to take (MCMC only)
#' @param ewald perform Wald CI computation instead of the default profile likelihood computation. This is the the 'FAST BMD' method of Ewald et al (2021)
#' @param transform Transforms doses using \eqn{\log(dose+\sqrt{dose^2+1})}. Note: this is a log transform that has a derivative defined when dose =0.
#' @return a model object
#' 
#' @examples 
#' M2           <- matrix(0,nrow=5,ncol=4)
#' colnames(M2) <- c("Dose","Resp","N","StDev")
#' M2[,1] <- c(0,25,50,100,200)
#' M2[,2] <- c(6,5.2,2.4,1.1,0.75)
#' M2[,3] <- c(20,20,19,20,20)
#' M2[,4] <- c(1.2,1.1,0.81,0.74,0.66)
#' model = single_continuous_fit(M2[,1,drop=F], M2[,2:4], BMD_TYPE="sd", BMR=1, ewald = T,
#'                              distribution = "normal",fit_type="laplace",model_type = "hill")
#' 
#' @export
single_continuous_fit <- function(D,Y,model_type="hill", fit_type = "laplace",
                                   prior=NA, BMD_TYPE = "sd", 
                                   BMR = 0.1, point_p = 0.01, distribution = "normal-ncv",
                                   alpha = 0.05, samples = 25000,degree=2,
                                   burnin = 1000,ewald = FALSE,
                                   transform = FALSE){
    Y <- as.matrix(Y) 
    D <- as.matrix(D) 

    myD = Y; 
    sstat = F # set sufficient statistics to false if there is only one column
    if (ncol(Y) > 1){
        sstat=T
    }
    
    if (!is.na(prior[1])){
      #parse the prior
      if (class(prior) != "BMD_Bayes_continuous_model"){
        stop("Prior is not the correct form. Please use a Bayesian Continuous Prior Model.")
      }
      t_prior_result <- parse_prior(prior)
      distribution <- t_prior_result$distribution
      model_type   <- t_prior_result$model
      prior = t_prior_result$prior
      PR = t_prior_result$prior
    }else{
      dmodel = which(model_type==.continuous_models)
      
      if (identical(dmodel, integer(0))){
        stop('Please specify one of the following model types: \n
            "hill","exp-3","exp-5","power","FUNL","polynomial"')
      }

      PR    = bayesian_prior_continuous_default(model_type,distribution,degree)
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
                     PR$priors[nrow(PR$priors),2]   = log(mean(Y[1,])/mean(Y[,3])^2)
                }
           }else{ 
                  if (distribution == "normal"){
                    PR$priors[nrow(PR$priors),2]   = log(var(Y[,1]))
                  }else{
                    PR$priors[nrow(PR$priors),2]   = log(mean(Y[,1])/var(Y[,1]))
                  }
           }
      }
      t_prior_result = create_continuous_prior(PR,model_type,distribution,degree)
      PR = t_prior_result$prior
    }
    dmodel = which(model_type==.continuous_models)

    type_of_fit = which(fit_type == c('laplace','mle','mcmc'))
    if (identical(type_of_fit,integer(0))){
      stop("Please choose one of the following fit types: 'laplace','mle','mcmc.' ")
    }
    
    rt = which(BMD_TYPE==c('abs','sd','rel','hybrid'))
    
    if (rt == 4){
      rt = 6; 
    }
    dis_type = which(distribution  == c("normal","normal-ncv","lognormal"))
    

    
    if(identical(dis_type,integer(0))){
      stop('Please specify the distribution as one of the following:\n
            "normal","normal-ncv","lognormal"')
    }
    
    DATA <- cbind(D,Y); 
    if (ncol(DATA)==4){
      colnames(DATA) =  c("Dose","Resp","N","StDev")
    }else if (ncol(DATA) == 2){
      colnames(DATA) =  c("Dose","Resp")
    }else{
      stop("The data do not appear to be in the correct format.")
    }
    #permute the matrix to the internal C values
    # Hill = 6, Exp3 = 3, Exp5 = 5, Power = 8, 
    # FUNL = 10
    permuteMat = cbind(c(1,2,3,4,5,6),c(6,3,5,8,10,666))
    fitmodel = permuteMat[dmodel,2]
    if (fitmodel == 10 && dis_type == 3){
         stop("The FUNL model is currently not defined for Log-Normal distribution.")
    }
    if (fitmodel == 666 && dis_type == 3){
         stop("Polynomial models are currently not defined for the Log-Normal distribution.")
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
      PR = MLE_bounds_continuous(model_type,distribution,degree, is_increasing)
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
      
      rvals <- run_continuous_single_mcmc(fitmodel,model_data$SSTAT,model_data$X,
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
        temp_me = temp_me[!is.infinite(temp_me[,1]),]
        temp_me = temp_me[!is.na(temp_me[,1]),]
        temp_me = temp_me[!is.nan(temp_me[,1]),]
        
        if( nrow(temp_me) > 5){
          te <- splinefun(temp_me[,2],temp_me[,1],method="hyman")
          rvals$bmd[2:3]  <- c(te(alpha),te(1-alpha))
        }else{
          rvals$bmd[2:3] <- c(NA,NA)
        }
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
      rvals   <- run_continuous_single(fitmodel,model_data$SSTAT,model_data$X,
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
        temp_me = temp_me[!is.infinite(temp_me[,1]),]
        temp_me = temp_me[!is.na(temp_me[,1]),]
        temp_me = temp_me[!is.nan(temp_me[,1]),]
        
        if( nrow(temp_me) > 5){
          te <- splinefun(temp_me[,2],temp_me[,1],method="hyman")
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

print.BMDcont_fit_MCMC<-function(p){
  
  BMDtype <- c('Absolute Deviation','Standard Deviation','Relative Deviation','Hybrid')
  
  
  cat ("Benchmark Dose Estimates using MCMC. \n")
  cat (sprintf("Continuous %s BMD: BMRF-%1.2f\n",BMDtype[p$options[1]],p$options[2]))
  cat (sprintf("Model Type: %s\n",p$model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  m <- mean(p$BMD)
  x <- quantile(p$BMD,c(p$options[4],1-p$options[4]))
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",m,x[1],x[2],100*(1-2*p$options[4]),"% 2-sided Confidence Interval"))
}

print.BMDcont_fit_laplace<-function(p){
  BMDtype <- c('Absolute Deviation','Standard Deviation','Relative Deviation','Hybrid')
  
  cat ("Benchmark Dose Estimates using Laplace \n")
  cat ("approximation to the Posterior\n")
  cat (sprintf("Continuous %s BMD: BMRF-%1.2f\n",BMDtype[p$options[1]],p$options[2]))
  cat (sprintf("Model Type: %s\n",p$model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",p$bmd[1],p$bmd[2],p$bmd[3],100*(1-2*p$options[4]),"% 2-sided Confidence Interval"))
}

print.BMDcont_fit_mle<-function(p){
  BMDtype <- c('Absolute Deviation','Standard Deviation','Relative Deviation','Hybrid')
  
  cat ("Benchmark Dose Estimates using MLE \n")
  cat (sprintf("Continuous %s BMD: BMRF-%1.2f\n",BMDtype[p$options[1]],p$options[2]))
  
  cat (sprintf("Model Type: %s\n",p$model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",p$bmd[1],p$bmd[2],p$bmd[3],100*(1-2*p$options[4]),"% 2-sided Confidence Interval"))
}



