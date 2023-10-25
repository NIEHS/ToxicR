.cont_exp_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * (1 - exp(-A[2] * doses^A[4])))
  return(out)
}
.cont_LMS_f <- function(A, doses, decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * (1 - exp(-A[2] * doses - A[4] * doses^2)))
}
.cont_invexp_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * (exp(-A[2] * doses^(-A[4]))))
  return(out)
}
.cont_gamma_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1)*pgamma(doses^A[4], shape = A[5], scale = 1.0 / A[2]))
  return(out)
}
.cont_gamma_efsa_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1)*pgamma(doses, shape = A[4], scale = 1.0 / A[2]))
  return(out)
}
.cont_invgamma_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1)*(1-pgamma(doses^(-A[4]), shape = A[5], scale = 1.0 / A[2])))
  return(out)
}
.cont_hill_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * (1 - A[2]^A[4] / (A[2]^A[4] + doses^A[4]))) 
  return(out)
}
.cont_lognormal_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * pnorm(log(A[2]) + A[4]*log(doses)))
  return(out)
}
.cont_lomax_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * (1 - (A[2] / (A[2] + doses^A[4]))^A[5]))
  return(out)
}
.cont_invlomax_aerts_f <- function(A,doses,decrease=F){
  out <- A[1] * (1 + (A[3] - 1) * ((A[2] / (A[2] + doses^(-A[4])))^(A[5])))
  return(out)
}
.cont_logskew_aerts_f <- function(A,doses,decrease=F){
  #out <- A[1] * (1 + (A[3] - 1) * psn(log(A[2]) + A[4]*log(doses), alpha=A[5]))
  temp <- sapply(doses, function(x) .owenst_fn(log(A[2]) + A[4]*log(x), A[5]))
  cdf <- pnorm(log(A[2]) + A[4]*log(doses)) - 2.0*temp
  out <- A[1] * (1 + (A[3] - 1) * cdf)
  return(out)
}
.cont_invlogskew_aerts_f <- function(A,doses,decrease=F){
  #out <- A[1] * (1 + (A[3] - 1) * (1 - psn(log(A[2]) - A[4]*log(doses), alpha=A[5])))
  temp <- sapply(doses, function(x) .owenst_fn(log(A[2]) - A[4]*log(x), A[5]))
  cdf <- pnorm(log(A[2]) - A[4]*log(doses)) - 2.0*temp
  out <- A[1] * (1 + (A[3] - 1) * (1 - cdf))
  return(out)
}
.cont_logistic_aerts_f <- function(A,doses,decrease=F){
  out <- A[3] / (1 + exp(-A[1] - A[2] * doses^A[4]))
  return(out)
}
.cont_probit_aerts_f <- function(A,doses,decrease=F){
  out <- A[3] * pnorm(A[1] + A[2] * doses^A[4])
  return(out)
}

.cont_polynomial_f <- function(A,doses,decrease=F){

  B <- as.matrix(A,ncol=1)
  X <- matrix(1,nrow = length(doses),ncol=length(A))
  for (ii in 2:nrow(B)){
    X[,ii] = doses^(ii-1)
  }
  return(X%*%B)
}

# FUNL
.cont_FUNL_f <- function(A,doses,decrease=F){
     b <- A[1] + A[2]*exp(-exp(A[6])*(doses-A[5])^2)*(1/(1+exp(-(doses-A[3])/A[4])))
     return(b)
}

#dichotomous hill
.cont_hill_f <- function(parms,d,decrease=F){
  g  <- parms[1] 
  nu <- parms[2]
  k  <- parms[3];
  n  <- parms[4]; 
  rval <- g + nu*d^n/(k^n+d^n)
  return (rval)
}
#dichotomous log-logistic
.cont_exp_5_f <- function(parms,d,decrease=F){
  g <- parms[1]
  b <- parms[2];
  c <- parms[3];
  e <- parms[4]; 
  rval <- g*(exp(c)-(exp(c)-1.0)*(exp(-(b*d)^e)))
  return (rval)
}

#
.cont_exp_3_f <-function(parms,d,decrease = TRUE){
  if (decrease){
    f_sign = -1; 
  }else{
    f_sign = 1; 
  }
  g <- parms[1]
  b <- parms[2]
  e <- parms[3] 
  rval <- g*exp(f_sign*(b*d)^e)
  return (rval)
}

.cont_power_f <-function(parms,d,decrease=F){
  g <- parms[1]; 
  b <- parms[2];
  a <- parms[3]; 
  rval <- g + b*d^a
  return (rval)
}

.plot.BMDcont_fit_MCMC<-function(x,...){
  fit = x
  Dose <- NULL
  temp_args = list(...)
  if (!exists("qprob",temp_args)){
    qprob = 0.05
  }else{
    qprob = temp_args$qprob
  }
  
  isLogNormal = (grepl("Log-Normal",fit$full_model) == 1)
     
  IS_tranformed = fit$transformed
  density_col="blueviolet"
  credint_col="azure2"
  BMD_DENSITY = T
  if(is.null(IS_tranformed)){ #NULL for MA, which cannot transform
    IS_tranformed = FALSE
  }
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
  
  
  data_d = fit$data
  
  IS_SUFFICIENT = FALSE
  if (ncol(data_d) == 4 ){ #sufficient statistics
    IS_SUFFICIENT = TRUE
    mean <- data_d[,2,drop=F]
    se   <- data_d[,4,drop=F]/sqrt(fit$data[,3,drop=F])
    doses = data_d[,1,drop=F]
    uerror <- mean+2*se
    lerror <- mean-2*se
    dose = c(doses,doses)
    Response = c(uerror,lerror)
    lm_fit = lm(mean ~ doses,weights = 1/se*se)
  }else{
    Response <- data_d[,2,drop=F]
    doses = data_d[,1,drop=F]
    lm_fit = lm(Response~doses)
  }
  
  if (coefficients(lm_fit)[2] < 0){
    decrease = TRUE
  }else{
    decrease = FALSE
  }
  
  # Single Model 
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/500)
  if( IS_tranformed)
  {
    test_doses = asinh(test_doses)
  }
  Q <- switch(fit$model,
              "exp-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_exp_aerts_f, d=test_doses,decrease=decrease),
              "invexp-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_invexp_aerts_f, d=test_doses,decrease=decrease),
              "gamma-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_gamma_aerts_f, d=test_doses,decrease=decrease),
              "invgamma-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_invgamma_aerts_f, d=test_doses,decrease=decrease),
              "hill-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_hill_aerts_f, d=test_doses,decrease=decrease),
              "lomax-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_lomax_aerts_f, d=test_doses,decrease=decrease), 
              "invlomax-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_invlomax_aerts_f, d=test_doses,decrease=decrease), 
              "lognormal-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_lognormal_aerts_f, d=test_doses,decrease=decrease), 
              "logskew-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_logskew_aerts_f, d=test_doses,decrease=decrease), 
              "invlogskew-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_invlogskew_aerts_f, d=test_doses,decrease=decrease), 
              "logistic-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_logistic_aerts_f, d=test_doses,decrease=decrease), 
              "probit-aerts" = apply(fit$mcmc_result$PARM_samples,1,.cont_probit_aerts_f, d=test_doses,decrease=decrease),
              "LMS" = apply(fit$mcmc_result$PARM_samples,1,.cont_LMS_f, d=test_doses,decrease=decrease),
              "gamma-efsa" = apply(fit$mcmc_result$PARM_samples,1,.cont_gamma_efsa_f, d=test_doses,decrease=decrease))
  if(isLogNormal & !is.null(Q)){
    Q <- exp(Q + exp(fit$mcmc_result$PARM_samples[,ncol(fit$mcmc_result$PARM_samples)])/2)
  }
  if (fit$model=="FUNL"){
     Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_FUNL_f, d=test_doses,decrease=decrease)   
  }
  if (fit$model=="hill"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_hill_f, d=test_doses,decrease=decrease)
  }
  if (fit$model=="exp-3"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_exp_3_f, d=test_doses,decrease=decrease)
    if (isLogNormal){
       Q <- exp(log(Q)+
        exp(fit$mcmc_result$PARM_samples[,ncol(fit$mcmc_result$PARM_samples)])/2)
    }
  }
  if (fit$model=="exp-5"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_exp_5_f, d=test_doses,decrease=decrease)
    if (isLogNormal){
         Q <- exp(log(Q)+
                       exp(fit$mcmc_result$PARM_samples[,ncol(fit$mcmc_result$PARM_samples)])/2)
    }
  }
  if (fit$model=="power"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_power_f, d=test_doses,decrease=decrease)
  }
  if (fit$model=="polynomial"){
    if (length(grep(": normal-ncv", tolower(fit$full_model)))>0){
      degree = ncol(fit$mcmc_result$PARM_samples) - 2
    }else{
      degree = ncol(fit$mcmc_result$PARM_samples) - 1
    }
    Q <- apply(fit$mcmc_result$PARM_samples[,1:degree],1,.cont_polynomial_f, 
               d=test_doses,decrease=decrease)
    
  }
  if( IS_tranformed)
  {
    test_doses = sinh(test_doses)
  }
  
  Q <- t(Q)
  me <- apply(Q,2,quantile, probs = 0.5)
  lq <- apply(Q,2,quantile, probs = qprob)
  uq <- apply(Q,2,quantile, probs = 1-qprob)
  
  # Continuous case density? 
  temp_fit <- splinefun(test_doses,me)
  ma_mean = temp_fit
  # Geom_polygon ? etc..
  plot_gg <- ggplot() +xlim(-max(test_doses)*5,max(test_doses)*5)+
    geom_line(aes(x=test_doses,y=me),color="blue",linewidth=2)+
    labs(x="Dose", y="Response",title=paste(fit$full_model, "MCMC",sep=",  Fit Type: " ))+
    theme_minimal()

  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
   
    plot_gg <- plot_gg +
      geom_segment(aes(x=fit$bmd[2], y=ma_mean(fit$bmd[1]), xend=fit$bmd[3],
                       yend=ma_mean(fit$bmd[1])),color="darkslategrey",linewidth=1.2, alpha=0.9) +
      annotate( geom = "text", x = fit$bmd[2], y = ma_mean(fit$bmd[1]),
                label = "[", size = 10,color="darkslategrey", alpha=0.9)+
      annotate(geom = "text", x = fit$bmd[3], y = ma_mean(fit$bmd[1]),
               label = "]", size = 10,color="darkslategrey", alpha=0.9) +
      annotate(geom = "point", x = fit$bmd[1], y = ma_mean(fit$bmd[1]),
               size = 5, color="darkslategrey",shape=17, alpha=0.9)
   }
  
# Add density 
  if (BMD_DENSITY ==TRUE){
    temp = fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp = temp[!is.na(temp)]
    temp = temp[!is.infinite(temp)]
    # Dens =  density(temp,cut=c(max(test_doses)), n=512, from=0, to=max(test_doses))
    # Dens =  density(temp,cut=c(max(test_doses)),adjust =1.5, n=512, from=min(test_doses), to=max(test_doses),na.rm=TRUE)
    errorFun <-function(e) {
      y <- rep(1,512)
      x <- seq(min(test_doses),max(test_doses), (max(test_doses) - min(test_doses))/511)
      return(data.frame(x=x,y=y))
    }

    # Dens =  density(temp,cut=c(max(test_doses)),adjust =1.5, n=512, from=min(test_doses), to=max(test_doses))
    Dens = tryCatch(density(temp,cut=c(max(test_doses)),adjust =1.5, n=512, 
                            from=min(test_doses), to=max(test_doses)), 
                    error = errorFun)
    Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.6
    temp = which(Dens$x < max(test_doses))
    D1_y = Dens$y[temp]
    D1_x = Dens$x[temp]
    qm = min(Response)
    scale = (max(Response)-min(Response))/max(D1_y) *.40
    plot_gg<-plot_gg +
             geom_polygon(aes(x=c(0,D1_x,max(doses)),y=c(min(Response),
                        min(Response)+D1_y*scale,min(Response))), fill = "blueviolet", alpha=0.6)
  }

  width=3
  
  if (IS_SUFFICIENT){
    plot_gg<- plot_gg +
      geom_errorbar(aes(x=doses, ymin=lerror, ymax=uerror),color="black",size=0.8,width=width)+
      geom_point(aes(x=doses,y=mean),size=3, shape=21, fill="white")
    
  }else{
    data_in<-data.frame(cbind(doses,Response))
    plot_gg<-plot_gg +
      geom_point(data=data_in,aes(x=doses,y=Response))
  }

  plot_gg <-plot_gg +
            geom_polygon(aes(x=c(test_doses,test_doses[length(test_doses):1]),y=c(uq,lq[length(test_doses):1])),fill="blue",alpha=0.1)
  return(plot_gg + coord_cartesian(xlim=c(min(test_doses),max(test_doses)),expand=F))

  
}
  
# This part matches with single_continous_fit part- SL 06/02/21 
.plot.BMDcont_fit_maximized<-function(x,...){
  A = x
  temp_args = list(...)
  if (!exists("qprob",temp_args)){
    qprob = 0.05
  }else{
    qprob = temp_args$qprob
  }
  
  isLogNormal = (grepl("Log-Normal",A$full_model) == 1)
  IS_tranformed = A$transformed
  fit <-A
  density_col="blueviolet"
  credint_col="azure2"

  if(is.null(IS_tranformed)){ #NULL for MA, which cannot transform
    IS_tranformed = FALSE
  }
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
     
  data_d = A$data
  IS_SUFFICIENT = FALSE
  
  if (ncol(data_d) == 4 ){ #sufficient statistics
    IS_SUFFICIENT = TRUE
    mean <- data_d[,2,drop=F]
    se   <- data_d[,4,drop=F]/sqrt(fit$data[,3,drop=F])
    doses = data_d[,1,drop=F]
    uerror <- mean+2*se
    lerror <- mean-2*se
    dose = c(doses,doses)
    Response = c(uerror,lerror)
    lm_fit = lm(mean ~ doses,weights = 1/(se*se))
  }else{
    Response <- data_d[,2,drop=F]
    doses = data_d[,1,drop=F]
    lm_fit = lm(Response~doses)
  }
  
  if (coefficients(lm_fit)[2] < 0){
    decrease = TRUE
  }else{
    decrease = FALSE
  }
  
  # I fixed some logic of inputs in if/else statement- they used to be fit$data
  # SL : Should Plot's x axis be based on test_dose? 
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)-min(doses))/500)
  if (IS_tranformed)
  {
    test_doses <- asinh(test_doses)
  }
  #Pre defined function- lm_fit can be used for fitting parameters?
  me <- switch(fit$model,
              "exp-aerts" = .cont_exp_aerts_f(fit$parameters,test_doses),
              "invexp-aerts" = .cont_invexp_aerts_f(fit$parameters,test_doses),
              "gamma-aerts" = .cont_gamma_aerts_f(fit$parameters,test_doses),
              "invgamma-aerts" = .cont_invgamma_aerts_f(fit$parameters,test_doses),
              "hill-aerts" = .cont_hill_aerts_f(fit$parameters,test_doses),
              "lomax-aerts" = .cont_lomax_aerts_f(fit$parameters,test_doses), 
              "invlomax-aerts" = .cont_invlomax_aerts_f(fit$parameters,test_doses), 
              "lognormal-aerts" = .cont_lognormal_aerts_f(fit$parameters,test_doses), 
              "logskew-aerts" = .cont_logskew_aerts_f(fit$parameters,test_doses), 
              "invlogskew-aerts" = .cont_invlogskew_aerts_f(fit$parameters,test_doses), 
              "logistic-aerts" = .cont_logistic_aerts_f(fit$parameters,test_doses), 
              "probit-aerts" = .cont_probit_aerts_f(fit$parameters,test_doses),
              "LMS" = .cont_LMS_f(fit$parameters, test_doses),
              "gamma-efsa" = .cont_gamma_efsa_f(fit$parameters, test_doses))
  if(isLogNormal & !is.null(me)){
    me <- exp(me)
  }
  if (fit$model=="FUNL"){
     me <- .cont_FUNL_f(fit$parameters,test_doses)
  }  
  if (fit$model=="hill"){
    me <- .cont_hill_f(fit$parameters,test_doses)
  }
  if (fit$model=="exp-3"){
    me <- .cont_exp_3_f(fit$parameters,test_doses,decrease)
  }
  if (fit$model=="exp-5"){
    me <- .cont_exp_5_f(fit$parameters,test_doses)
  }
  if (fit$model=="power"){
    me <- .cont_power_f(fit$parameters,test_doses)
  }
  if (fit$model=="polynomial"){
    if (length(grep(": normal-ncv", tolower(fit$full_model)))>0){
      degree = length(fit$parameters) - 2
    }else{
      degree = length(fit$parameters) - 1
    }
   
    me <- .cont_polynomial_f(fit$parameters[1:degree],test_doses)
  }
  if (isLogNormal){
       var = exp(fit$parameters[length(fit$parameters)])
       me = exp(log(me)+var/2)
  }
  
  if (IS_tranformed)
  {
    test_doses <- sinh(test_doses)
  }
  
  temp_fit <- splinefun(test_doses,me)
  ma_mean  <- temp_fit
  plot_gg<-ggplot()+
          geom_line(aes(x=test_doses,y=me),color="blue",linewidth=2)+xlim(-max(test_doses)*5,max(test_doses)*5)+
          labs(x="Dose", y="Response",title=paste(fit$full_model, "Maximized",sep=",  Fit Type: " ))+
          theme_minimal()
        

  
  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    if (!sum(is.na(fit$bmd))){
        
      plot_gg <- plot_gg +
        geom_segment(aes(x=fit$bmd[2], y=ma_mean(fit$bmd[1]), xend=fit$bmd[3],
                         yend=ma_mean(fit$bmd[1])),color="darkslategrey",linewidth=1.2, alpha=0.9) +
        annotate( geom = "text", x = fit$bmd[2], y = ma_mean(fit$bmd[1]),
                  label = "[", size = 10,color="darkslategrey", alpha=0.9)+
        annotate(geom = "text", x = fit$bmd[3], y = ma_mean(fit$bmd[1]),
                 label = "]", size = 10,color="darkslategrey", alpha=0.9) +
        annotate(geom = "point", x = fit$bmd[1], y = ma_mean(fit$bmd[1]),
                 size = 5, color="darkslategrey",shape=17, alpha=0.9)
    }
  }
  # Assign them temporarily 
  width=3
  
  if (IS_SUFFICIENT){
    plot_gg<- plot_gg +
              geom_errorbar(aes(x=doses, ymin=lerror, ymax=uerror),color="grey",size=0.5, width=3)+
              geom_point(aes(x=doses,y=mean),size=3, shape=21, fill="white")
        
  }else{
    data_in<-data.frame(cbind(doses,Response))
    plot_gg<-plot_gg +
          geom_point(aes(x=doses,y=Response))
  }
  
  return(plot_gg + coord_cartesian(xlim=c(min(test_doses),max(test_doses)),expand=F))
  
  
}

# Base plot- MCMC or BMD?
.plot.BMDcontinuous_MA <- function(x,...){
  A = x
  model_no <- x_axis <- y_axis <-cols <- NULL
  temp_args = list(...)
  if (!exists("qprob",temp_args)){
    qprob = 0.05
  }else{
    qprob = temp_args$qprob
  }
  # Should be matched with BMD_MA plots
  # SL 06/02 Updated 
  # Later, we'll have it 
     density_col="blueviolet"
     credint_col="azure2"
     class_list <- names(A)
     fit_idx    <- grep("Indiv_",class_list)
  
     A$posterior_probs[!is.finite( A$posterior_probs)] = 0.0
     #plot the model average curve
     if ("BMDcontinuous_MA_mcmc" %in% class(A)){ # mcmc run
          n_samps <- nrow(A[[fit_idx[1]]]$mcmc_result$PARM_samples)
          data_d   <-  A[[fit_idx[1]]]$data
          max_dose <- max(data_d[,1])
          min_dose <- min(data_d[,1])
          test_doses <- seq(min_dose,max_dose,(max_dose-min_dose)/500) 
          ma_samps <- sample(fit_idx,n_samps, replace=TRUE,prob = A$posterior_probs)
          temp_f   <- matrix(0,n_samps,length(test_doses))
          temp_bmd <- rep(0,length(test_doses))
          width= (max_dose-min_dose)/20
          
          # 06/07/21 SL Update
          IS_SUFFICIENT=FALSE
   
          if (ncol(data_d) == 4 ){ #sufficient statistics
            IS_SUFFICIENT = TRUE
            mean <- data_d[,2,drop=F]
            se   <- data_d[,4,drop=F]/sqrt(data_d[,3,drop=F])
            doses = data_d[,1,drop=F]
            uerror <- mean+2*se
            lerror <- mean-2*se
            dose = c(doses,doses)
            Response = c(uerror,lerror)
            lm_fit = lm(mean ~ doses,weights = 1/(se*se))
          }else{
            Response <- data_d[,2,drop=F]
            doses = data_d[,1,drop=F]
            lm_fit = lm(Response~doses)
          }
          
          if (coefficients(lm_fit)[2] < 0){
            decrease = TRUE
          }else{
            decrease = FALSE
          }
          
          for (ii in 1:n_samps){
               fit <- A[[fit_idx[ma_samps[ii]]]]
               isLogNormal = (grepl("Log-Normal",fit$full_model) == 1)
               asdf <- switch(fit$model,
                           "exp-aerts" = .cont_exp_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "invexp-aerts" = .cont_invexp_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "gamma-aerts" = .cont_gamma_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "invgamma-aerts" = .cont_invgamma_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "hill-aerts" = .cont_hill_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "lomax-aerts" = .cont_lomax_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses), 
                           "invlomax-aerts" = .cont_invlomax_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses), 
                           "lognormal-aerts" = .cont_lognormal_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses), 
                           "logskew-aerts" = .cont_logskew_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses), 
                           "invlogskew-aerts" = .cont_invlogskew_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses), 
                           "logistic-aerts" = .cont_logistic_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses), 
                           "probit-aerts" = .cont_probit_aerts_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "LMS" = .cont_LMS_f(fit$mcmc_result$PARM_samples[ii,],test_doses),
                           "gamma-efsa" = .cont_gamma_efsa_f(fit$mcmc_result$PARM_samples[ii,],test_doses))
               if(isLogNormal & !is.null(asdf)){
                 asdf <- exp(asdf)
               }
               if(!is.null(asdf)){
                 temp_f[ii,] <- asdf
                 temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
                 rm(asdf)
               }
               if (fit$model=="FUNL"){
                    temp_f[ii,] <- .cont_FUNL_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }  
               if (fit$model=="hill"){
                    temp_f[ii,] <- .cont_hill_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
               if (fit$model=="exp-3"){
                    temp_f[ii,] <- .cont_exp_3_f(fit$mcmc_result$PARM_samples[ii,],test_doses,decrease)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
               if (fit$model=="exp-5"){
                    temp_f[ii,] <- .cont_exp_5_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
               if (fit$model=="power"){
                    temp_f[ii,] <- .cont_power_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
          }
          
          temp_f[is.infinite(temp_f)] = NA
          temp_f[abs(temp_f) > 1e10] = NA
          
          # If temp_bmd== Inf then delete;
          # Updated 06/02/21 SL
          temp_bmd[is.infinite(temp_bmd)] = NA
          
          me <- apply(temp_f,2,quantile, probs = 0.5,na.rm = TRUE) # BMD
          lq <- apply(temp_f,2,quantile, probs = qprob,na.rm = TRUE) # BMDL
          uq <- apply(temp_f,2,quantile, probs = 1-qprob,na.rm = TRUE) # BMDU

          
          # 06/02/21 SL update
          if (IS_SUFFICIENT){
         
              plot_gg<-ggplot()+xlim(-max(test_doses)*5,min(test_doses)*5)+
                  geom_point(aes(x=data_d[,1],y=data_d[,2]))+
                  geom_errorbar(aes(x=data_d[,1], ymin=lerror, ymax=uerror),color="grey",size=0.8,width=width)+
                  xlim(c(min(data_d[,1])-width,max(data_d[,1])*1.03))+
                  labs(x="Dose", y="Response",title="Continous MA fitting")+
                  theme_minimal()
          
          }else{
            plot_gg<-ggplot()+xlim(-max(test_doses)*5,min(test_doses)*5)+
              geom_point(aes(x=doses,y=Response))+
              xlim(c(min(doses),max(doses)*1.03))+
              labs(x="Dose", y="Response",title="Continous MA fitting")+
              theme_minimal()
          }
          
          
          plot_gg<-plot_gg+
                   geom_ribbon(aes(x=test_doses,ymin=lq,ymax=uq),fill="blue",alpha=0.1)
          
          plot_gg<-plot_gg+
                   geom_line(aes(x=test_doses,y=me),col="blue",linewidth=2)
         
          bmd <- quantile(temp_bmd,c(qprob,0.5,1-qprob),na.rm = TRUE)
  
          
          
          if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
            temp = temp_bmd[temp_bmd < 10*max(test_doses)]
            temp = temp[!is.infinite(temp_bmd)]
            temp = temp[!is.na(temp)]
       
            # Density only creates few data points SL
            
            # Fixed part 06/04/21
            Dens =  density(temp,cut=c(5*max(test_doses)), n=1000, from=0, to=max(test_doses))
          
            Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.6
            temp = which(Dens$x < max(test_doses))
            D1_y = Dens$y[temp]
            D1_x = Dens$x[temp]
            qm = min(Response)
            scale = (max(Response)-min(Response))/max(D1_y) *.40
            
          
             plot_gg<-plot_gg+
                    geom_polygon(aes(x=c(max(0,min(D1_x)),D1_x,max(D1_x)),
                                     y=c(min(Response),min(Response)+D1_y*scale,min(Response))),
                                     fill = "blueviolet", alpha=0.6)

          }
          
          ## 
          # Add lines to the BMD
          ma_mean <- splinefun(test_doses,me)
          ma_BMD = A$bmd
       
          df<-NULL
           
          # Problem of the loop using this case- the ggplot is not added automatically, 
          # It replaces the last one;
          
          for (ii in 1:length(fit_idx)){
            if (!is.finite(A$posterior_probs[ii])){
              A$posterior_probs[ii] = 0
            }
            if (A$posterior_probs[ii]>0.05){
               fit <- A[[fit_idx[ii]]]
               isLogNormal = (grepl("Log-Normal",fit$full_model) == 1)
               f <- switch(fit$model,
                            "exp-aerts" = .cont_exp_aerts_f(fit$parameters,test_doses),
                            "invexp-aerts" = .cont_invexp_aerts_f(fit$parameters,test_doses),
                            "gamma-aerts" = .cont_gamma_aerts_f(fit$parameters,test_doses),
                            "invgamma-aerts" = .cont_invgamma_aerts_f(fit$parameters,test_doses),
                            "hill-aerts" = .cont_hill_aerts_f(fit$parameters,test_doses),
                            "lomax-aerts" = .cont_lomax_aerts_f(fit$parameters,test_doses), 
                            "invlomax-aerts" = .cont_invlomax_aerts_f(fit$parameters,test_doses), 
                            "lognormal-aerts" = .cont_lognormal_aerts_f(fit$parameters,test_doses), 
                            "logskew-aerts" = .cont_logskew_aerts_f(fit$parameters,test_doses), 
                            "invlogskew-aerts" = .cont_invlogskew_aerts_f(fit$parameters,test_doses), 
                            "logistic-aerts" = .cont_logistic_aerts_f(fit$parameters,test_doses), 
                            "probit-aerts" = .cont_probit_aerts_f(fit$parameters,test_doses),
                            "LMS" = .cont_LMS_f(fit$parameters,test_doses),
                            "gamma-efsa" = .cont_gamma_efsa_f(fit$parameters,test_doses))
               if(isLogNormal & !is.null(f)){
                 f <- exp(f)
               }
               if (fit$model=="FUNL"){
                    f <- .cont_FUNL_f(fit$parameters,test_doses)
               }  
               if (fit$model=="hill"){
                    f <- .cont_hill_f(fit$parameters,test_doses)
               }
               if (fit$model=="exp-3"){
                   temp = fit$parameters
                    f <- .cont_exp_3_f(temp,test_doses,decrease)
               }
               if (fit$model=="exp-5"){
                    f <- .cont_exp_5_f(fit$parameters,test_doses)
               }
               if (fit$model=="power"){
                    f <- .cont_power_f(fit$parameters,test_doses)
               }
               col = 'coral3'
               temp_df<-data.frame(x_axis=test_doses,y_axis=f,cols=col,model_no=ii, alpha_lev=unname(A$posterior_probs[ii]))
               # # 06/19/21 SL update 
               df     <-data.frame(x_axis=test_doses,y_axis=f,cols=col,model_no=ii, alpha_lev=unname(A$posterior_probs[ii]))
               
               df <-rbind(df,temp_df)
        
               #SL Updated 06/18/21 -- Transparency update based on posterior probability and Y scale for dichotomous case
               temp_data<- df %>% filter(model_no==ii)
               
               plot_gg<- plot_gg+
                 geom_line(data=temp_data, aes(x=x_axis,y=y_axis,color=cols),alpha=unique(temp_data$alpha_lev),show.legend=F)+
                 theme_minimal()
               
               plot_gg <- plot_gg +
                         geom_segment(aes(x=A$bmd[2], y=ma_mean(A$bmd[1]), xend=min(max(doses),A$bmd[3]),
                                          yend=ma_mean(A$bmd[1])),color="darkslategrey",linewidth=1.2, alpha=0.9) +
                         annotate( geom = "text", x = A$bmd[2], y = ma_mean(A$bmd[1]),
                                   label = "[", size = 10,color="darkslategrey", alpha=0.9)+
                         annotate(geom = "text", x = A$bmd[3], y = ma_mean(A$bmd[1]),
                                  label = "]", size = 10,color="darkslategrey", alpha=0.9) +
                         annotate(geom = "point", x = A$bmd[1], y = ma_mean(A$bmd[1]),
                                  size = 5, color="darkslategrey",shape=17, alpha=0.9)
            }
            
          
          }
          
          

     }else{ #laplace run
       
       data_d   <-  A[[fit_idx[1]]]$data
       max_dose <- max(data_d[,1])
       min_dose <- min(data_d[,1])
       width= (max_dose-min_dose)/20
       test_doses <- seq(min_dose,max_dose,(max_dose-min_dose)/200); 
       temp_bmd <- rep(0,length(test_doses))
       IS_SUFFICIENT = F
       if (ncol(data_d) == 4 ){ #sufficient statistics
         mean <- data_d[,2,drop=F]
         se   <- data_d[,4,drop=F]/sqrt(data_d[,3,drop=F])
         doses = data_d[,1,drop=F]
         uerror <- mean+2*se
         lerror <- mean-2*se
         dose = c(doses,doses)
         Response = c(uerror,lerror)
       
         lm_fit = lm(mean ~ doses,weights = 1/(se*se))
         IS_SUFFICIENT = T
       }else{
         Response <- data_d[,2,drop=F]
         doses = data_d[,1,drop=F]
         lm_fit = lm(Response~doses)
       }
       
       
       if (coefficients(lm_fit)[2] < 0){
         decrease = TRUE
       }else{
         decrease = FALSE
       }
       me = test_doses*0   
       for (ii in 1:length(fit_idx)){
         fit <- A[[fit_idx[ii]]]
         isLogNormal = (grepl("Log-Normal",fit$full_model) == 1)
         t <- switch(fit$model,
                      "exp-aerts" = .cont_exp_aerts_f(fit$parameters,test_doses),
                      "invexp-aerts" = .cont_invexp_aerts_f(fit$parameters,test_doses),
                      "gamma-aerts" = .cont_gamma_aerts_f(fit$parameters,test_doses),
                      "invgamma-aerts" = .cont_invgamma_aerts_f(fit$parameters,test_doses),
                      "hill-aerts" = .cont_hill_aerts_f(fit$parameters,test_doses),
                      "lomax-aerts" = .cont_lomax_aerts_f(fit$parameters,test_doses), 
                      "invlomax-aerts" = .cont_invlomax_aerts_f(fit$parameters,test_doses), 
                      "lognormal-aerts" = .cont_lognormal_aerts_f(fit$parameters,test_doses), 
                      "logskew-aerts" = .cont_logskew_aerts_f(fit$parameters,test_doses), 
                      "invlogskew-aerts" = .cont_invlogskew_aerts_f(fit$parameters,test_doses), 
                      "logistic-aerts" = .cont_logistic_aerts_f(fit$parameters,test_doses), 
                      "probit-aerts" = .cont_probit_aerts_f(fit$parameters,test_doses),
                      "LMS" = .cont_LMS_f(fit$parameters,test_doses),
                      "gamma-efsa" = .cont_gamma_efsa_f(fit$parameters,test_doses))
         if (!is.finite(A$posterior_probs[ii])){
           A$posterior_probs[ii] = 0
         }
         if(isLogNormal & !is.null(t)){
           t <- exp(t)
         }
         if(!is.null(t)){
           if (!is.finite(A$posterior_probs[ii])){
             A$posterior_probs[ii] = 0
           }
           if(A$posterior_probs[ii] > 0){
             me = t*A$posterior_probs[ii] + me
           }
         }
         if (fit$model=="FUNL"){
           t <- .cont_FUNL_f(fit$parameters,test_doses)
           if (!is.finite(A$posterior_probs[ii])){
             A$posterior_probs[ii] = 0
           }
           if(A$posterior_probs[ii] > 0){
             me = t*A$posterior_probs[ii] + me
           }
          
         }  
         if (fit$model=="hill"){
            
           t <- .cont_hill_f(fit$parameters,test_doses)
           
           # SL comment - why the name of object is BB? At the beginning it was declared as A-  05/28/21
           # I guess this part should be A as well 
           if (!is.finite(A$posterior_probs[ii])){
             A$posterior_probs[ii] = 0
           }
           if(A$posterior_probs[ii] > 0){
             me = t*A$posterior_probs[ii] + me
           }
         }
         if (fit$model=="exp-3"){
           t <- .cont_exp_3_f(fit$parameters,test_doses,decrease)
          
           if(A$posterior_probs[ii] > 0){
             me = t*A$posterior_probs[ii] + me
           }
         }
         if (fit$model=="exp-5"){
           t <- .cont_exp_5_f(fit$parameters,test_doses)
           if(A$posterior_probs[ii] > 0){
             me = t*A$posterior_probs[ii] + me
           }
         }
         if (fit$model=="power"){
           t <- .cont_power_f(fit$parameters,test_doses)
           if(A$posterior_probs[ii] > 0){
             me = t*A$posterior_probs[ii] + me
           }
         }
       }

       if (IS_SUFFICIENT){

         plot_gg<-ggplot()+xlim(-max(test_doses)*5,min(test_doses)*5)+
           geom_point(aes(x=data_d[,1],y=data_d[,2]))+
           geom_errorbar(aes(x=data_d[,1], ymin=lerror, ymax=uerror),color="grey",size=0.8,width=width)+
           xlim(c(min(data_d[,1])-width,max(data_d[,1])*1.03))+
           labs(x="Dose", y="Response",title="Continous MA fitting")+
           theme_minimal()
           y_min = min(lerror)
           y_max = max(uerror)
       }else{
         plot_gg<-ggplot()+xlim(-max(test_doses)*5,min(test_doses)*5)+
           geom_point(aes(x=doses,y=Response))+
           xlim(c(min(doses),max(doses)*1.03))+
           labs(x="Dose", y="Response",title="Continous MA fitting")+
           theme_minimal()
           y_min = min(Response)
           y_max = max(Response)
       }
       
        
       plot_gg<-plot_gg+
         geom_line(aes(x=test_doses,y=me),col="blue",linewidth=2)
 
       ## 
       # Add lines to the BMD
       ma_mean <- splinefun(test_doses,me)
       ma_BMD = A$bmd

       # Not sure about this part - SL 05/28/21
       #Plot only level >2
       
       
       for (ii in 1:length(fit_idx)){
       df<-NULL
         if (!is.finite(A$posterior_probs[ii])){
           A$posterior_probs[ii] = 0
         }
         if (A$posterior_probs[ii]>0.05){
           fit <- A[[fit_idx[ii]]]
           isLogNormal = (grepl("Log-Normal",fit$full_model) == 1)
           f <- switch(fit$model,
                       "exp-aerts" = .cont_exp_aerts_f(fit$parameters,test_doses),
                       "invexp-aerts" = .cont_invexp_aerts_f(fit$parameters,test_doses),
                       "gamma-aerts" = .cont_gamma_aerts_f(fit$parameters,test_doses),
                       "invgamma-aerts" = .cont_invgamma_aerts_f(fit$parameters,test_doses),
                       "hill-aerts" = .cont_hill_aerts_f(fit$parameters,test_doses),
                       "lomax-aerts" = .cont_lomax_aerts_f(fit$parameters,test_doses), 
                       "invlomax-aerts" = .cont_invlomax_aerts_f(fit$parameters,test_doses), 
                       "lognormal-aerts" = .cont_lognormal_aerts_f(fit$parameters,test_doses), 
                       "logskew-aerts" = .cont_logskew_aerts_f(fit$parameters,test_doses), 
                       "invlogskew-aerts" = .cont_invlogskew_aerts_f(fit$parameters,test_doses), 
                       "logistic-aerts" = .cont_logistic_aerts_f(fit$parameters,test_doses), 
                       "probit-aerts" = .cont_probit_aerts_f(fit$parameters,test_doses),
                       "LMS" = .cont_LMS_f(fit$parameters,test_doses),
                       "gamma-efsa" = .cont_gamma_efsa_f(fit$parameters,test_doses))
           if(isLogNormal & !is.null(f)){
             f <- exp(f)
           }
           if (fit$model=="FUNL"){
             f <- .cont_FUNL_f(fit$parameters,test_doses)
           }
           if (fit$model=="hill"){
             f <- .cont_hill_f(fit$parameters,test_doses)
           }
           if (fit$model=="exp-3"){
             temp = fit$parameters
             f <- .cont_exp_3_f(temp,test_doses,decrease)
           }
           if (fit$model=="exp-5"){
             f <- .cont_exp_5_f(fit$parameters,test_doses)
           }
           if (fit$model=="power"){
             f <- .cont_power_f(fit$parameters,test_doses)
           }

           col = 'coral3'
           # Not using loop, but save data in the external data and load it later
           temp_df<-data.frame(x_axis=test_doses,y_axis=f,cols=col,model_no=ii, alpha_lev=unname(A$posterior_probs[ii]))
           df<-temp_df #rbind(df,temp_df)
           plot_gg<- plot_gg+
                geom_line(data=df, aes(x=x_axis,y=y_axis,color=col),alpha=0.5,show.legend=F)
         }
       

         
                    
       }
       plot_gg <- plot_gg +
                   geom_segment(aes(x=A$bmd[2], y=ma_mean(A$bmd[1]), xend=min(max(doses),abs(A$bmd[3])),
                                    yend=ma_mean(A$bmd[1])),color="darkslategrey",linewidth=1.2, alpha=0.9) +
                   annotate( geom = "text", x = A$bmd[2], y = ma_mean(A$bmd[1]),
                             label = "[", size = 10,color="darkslategrey", alpha=0.9)+
                   annotate(geom = "text", x = A$bmd[3], y = ma_mean(A$bmd[1]),
                            label = "]", size = 10,color="darkslategrey", alpha=0.9) +
                   annotate(geom = "point", x = A$bmd[1], y = ma_mean(A$bmd[1]),
                            size = 5, color="darkslategrey",shape=17, alpha=0.9)
     }
   
     return(plot_gg + 
              coord_cartesian(xlim=c(min(test_doses)-width,max(test_doses)+width),expand=F))
}
 