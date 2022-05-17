
#' Create a density plot from a model averaged model fit with MCMC.
#'
#' @title MAdensity_plot - Create a density plot from a model averaged model.
#' @param A the model averaged model to plot
#' @return Returns a \code{ggplot2} graphics object. 
#' @examples 
#' \donttest{
#' doses  <- cbind(c(0,25,50,100,200))
#' y <- cbind(c(6,5.2,2.4,1.1,0.75),
#'                 c(20,20,19,20,20),
#'                 c(1.2,1.1,0.81,0.74,0.66))
#' model <- ma_continuous_fit(doses,y,
#'                         fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1)
#' MAdensity_plot(model)
#' }
#' @export
MAdensity_plot <- function (A){
  #source("dicho_functions.R")
  UseMethod("MAdensity_plot")
}

# Sample Dichotomous Data set
.plot.density.BMDdichotomous_MA_MCMC<-function(A){
# Construct bmd sample plots for mcmc
  X1 <- X2 <- X3 <- NULL
  class_list <- names(A)
  fit_idx    <- grep("Individual_Model",class_list)
  qprob=0.05
  
  #Dose levels
  data<-A$Individual_Model_1$data
  doses<-data[,1]
  
  
  
  t_combine<-NA
  
  
  for (i in fit_idx){
    # Loop for the model
    fit<-A[[i]]
    test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
    probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])

    
    
    if (fit$model=="hill"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_hill_f, d=test_doses)
      
    }
    if (fit$model=="gamma"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_gamma_f, d=test_doses)
      
    }
    if (fit$model=="logistic"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_logist_f, d=test_doses)
      
    }
    if (fit$model=="log-logistic"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_llogist_f, d=test_doses)
      
    }
    if (fit$model=="probit"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_probit_f, d=test_doses)
      
    }
    if (fit$model=="log-probit"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_lprobit_f, d=test_doses)
      
    }
    if (fit$model=="multistage"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_multistage_f, d=test_doses)
      
    }
    if (fit$model=="qlinear"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_qlinear_f, d=test_doses)
    }
    
    temp <- fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp <- temp[!is.infinite(temp)]
    
    
    
    # Q <- t(Q)
    # 
    # me <- colMeans(Q)
    # lq <- apply(Q,2,quantile, probs = qprob)
    # uq <- apply(Q,2,quantile, probs = 1-qprob)
    
    Dens =  density(temp,cut=c(max(doses)))
    # what is this 0.4 means? Scale?
    
    # normalize it?-- We don't need it actually here
    # Dens$y = Dens$y/max(Dens$y) * max(probs)
    # temp = which(Dens$x < max(doses))
    # D1_y = Dens$y[temp]
    # D1_x = Dens$x[temp]
    
    
    
    # Do I need to stack up the dataset?
    
    
    temp_density<-data.frame(matrix(0,length(temp),3))
    temp_density[,2]=fit$model
    temp_density[,1]=temp
    temp_density[,3]=A$posterior_probs[i]
    
    # assign(paste("t",i,sep="_"),temp_density)
    # 06/21/21 Update
    t<-temp_density
    
    t_combine<-rbind(t_combine,t)
  }
  
  t_combine<-t_combine[-1,]
  
  # 
  # t_combine<-rbind(t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9)
  # 
  
  
  # This part is needed to get MA density plots
  # 
  # idx <- sample(1:9, length(A$Individual_Model_1$mcmc_result$BMD_samples),replace=TRUE,prob=A$posterior_probs)

  idx <- sample(1:length(fit_idx), length(A$Individual_Model_1$mcmc_result$BMD_samples),replace=TRUE,prob=A$posterior_probs)
  
  df<-NA
  # 
  # m1<-A$Individual_Model_1
  # c1<-m1$mcmc_result$BMD_samples
  # 
  # 
  # m2<-A$Individual_Model_2
  # c2<-m2$mcmc_result$BMD_samples
  # 
  # m3<-A$Individual_Model_3
  # c3<-m3$mcmc_result$BMD_samples
  # 
  # m4<-A$Individual_Model_4
  # c4<-m4$mcmc_result$BMD_samples
  # 
  # m5<-A$Individual_Model_5
  # c5<-m5$mcmc_result$BMD_samples
  # 
  # m6<-A$Individual_Model_6
  # c6<-m6$mcmc_result$BMD_samples
  # 
  # m7<-A$Individual_Model_7
  # c7<-m7$mcmc_result$BMD_samples
  # 
  # m8<-A$Individual_Model_8
  # c8<-m5$mcmc_result$BMD_samples
  # 
  # m9<-A$Individual_Model_9
  # c9<-m9$mcmc_result$BMD_samples
  # 
  # combine_samples<-data.frame(cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9))

  
  # How should I initialize this?
  for (i in 1:length(fit_idx)){
    m<-A[[i]]
    c<-m$mcmc_result$BMD_samples
    df<-data.frame(cbind(df,c))
  }
  
  df_samples<-data.frame(df[,-1])
  
  
  
  
  
  
  # Select MA values
  BMD_MA<-matrix(NA,length(A$Individual_Model_1$mcmc_result$BMD_samples),1)
  
  for (i in 1:length(A$Individual_Model_1$mcmc_result$BMD_samples)){
    # BMD_MA[i,1]<-combine_samples[sample(nrow(combine_samples), size=1, replace=TRUE),idx[i]]
    j<-sample(nrow(df_samples), size=1, replace=TRUE)
    BMD_MA[i,1]<-df_samples[j,idx[i]]  
    
  }
  
  BMD_MA<-data.frame(BMD_MA)
  
  t_ma<-BMD_MA %>%
    mutate(X1=BMD_MA,X2="Model Average",X3=1)
  #BMD_CDF should be shown here - it 
  t_ma2<-t_ma %>%
    select(X1,X2,X3)
  
  t_combine2<-rbind(t_combine,t_ma2)
  t_combine3<-t_combine2 %>%
    filter(as.numeric(X3)>0.05  , as.numeric(X1)!=Inf)
  
  
  
  
  # # Update to change color
  # 
  # step1<- p<-ggplot()+
  #   stat_density_ridges(data=t_combine3, aes(x=X1, y=fct_reorder(X2,X3,.desc=T),group=X2,
  #                                            fill=cut(X3,c(0,0.999,1))), 
  #                       calc_ecdf=TRUE,quantiles=c(0.025,0.975),na.rm=T,quantile_lines=T,scale=0.9)+
  #                       scale_fill_manual(name="X3", values=c("(0,0.999]"="grey", "(0.999,1]"="blue"))
  #   
  # 
  # step2<-step1+
  #   xlim(c(0,quantile(t_combine$X1,0.999)))+
  #   geom_vline(xintercept = A$bmd[1],linetype="longdash")+
  #   labs(x="Dose Level (Dotted Line : MA BMD)",y="", title="Density plots for each fitted model (Fit type: MCMC)")+theme_classic()
  # 
  # step2
  # 
  # step3<-step2+
  #   stat_density_ridges(data=t_combine3, aes(x=X1, y=fct_reorder(X2,X3,.desc=T), group=X2, fill = factor(stat(quantile))),
  #                                  geom = "density_ridges_gradient",
  #                                  calc_ecdf = TRUE, quantiles = c(0.025, 0.975), scale=0.9)+
  #   scale_fill_manual(
  #     name = "Probability",
  #     values = c("NA","NA","NA", "NA", "NA")
  #    )
  #   
  # 
  # step3
  # 

  
  p<-ggplot()+
    stat_density_ridges(data=t_combine3, aes(x=X1, y=fct_reorder(X2,X3,.desc=T),group=X2 ,alpha=sqrt(X3), fill=cut(X3,c(0,0.99,1))), 
                        calc_ecdf=TRUE,quantiles=c(0.025,0.975),na.rm=T,quantile_lines=T,scale=0.9)+
    xlim(c(0,quantile(t_combine$X1,0.99)))+
    geom_vline(xintercept = A$bmd[1],linetype="longdash")+
    scale_fill_manual(name="X3", values=c("(0,0.99]"="darkgrey", "(0.99,1]"="red"))+
    labs(x="Dose Level (Dotted Line : MA BMD)",y="", title="Density plots for each fitted model (Fit type: MCMC)")+
    theme_classic()
  
  

  p2<-p+theme(legend.position="none", axis.text.y=element_text(size=12))
    # stat_density_ridges(data=t_combine3, aes(x=X1, y=fct_reorder(X2,X3,.desc=T), group=X2,fill = factor(stat(quantile))),
    #                     geom = "density_ridges_gradient",
    #                     calc_ecdf = TRUE, quantiles = c(0.025, 0.975),scale=0.9
    # )+ scale_fill_manual(
    #   name = "factor(stat(quantile))",
    #   values = c("#FF0000A0", "NA", "#FF0000A0"),
    #   labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
    # )+
  
  
  p2
  return(p2) # Return output
  
  
}


# Laplace case- Is it even possible to do this? 
# No we don't need this part

.plot.density.BMDdichotomous_MA_maximized<-function(A){
  t_1 <- t_2 <- t_3 <- t_4 <- t_5 <- t_6 <- t_7 <- t_8 <- t_9 <- c3 <- X1 <- X2 <- X3 <- NULL
  class_list <- names(A)
  
  if (class(A)[2]=="BMDdichotomous_MA_maximized"){
    fit_idx <- grep("Fitted_Model",class_list)
    qprob=0.05
    
    #Dose levels
    data<-A$Fitted_Model_1$data
    doses<-data[,1]
    
  }else{
  fit_idx <- grep("Individual_Model",class_list)
  qprob=0.05
  
  #Dose levels
  data<-A$Individual_Model_1$data
  doses<-data[,1]
  }
  
  for (i in fit_idx){
    fit<-A[[i]]
    test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
    
    if (fit$model=="hill"){
      me <- .dich_hill_f(fit$parameters, d=test_doses)
    }else if (fit$model=="gamma"){
      me <- .dich_gamma_f(fit$parameters, d=test_doses)
    }else if (fit$model=="logistic"){
      me <- .dich_logist_f(fit$parameters, d=test_doses)
    }else if (fit$model=="log-logistic"){
      me <- .dich_llogist_f(fit$parameters, d=test_doses)    
    }else if (fit$model=="probit"){
      me <- .dich_probit_f(fit$parameters, d=test_doses)    
    }else if (fit$model=="log-probit"){
      me <- .dich_lprobit_f(fit$parameters, d=test_doses)    
    }else if (fit$model=="multistage"){
      me <- .dich_multistage_f(fit$parameters, d=test_doses)    
    }else if (fit$model=="qlinear"){
      me <- .dich_qlinear_f(fit$parameters, d=test_doses)    
    }else if (fit$model=="weibull"){
      me <- .dich_weibull_f(fit$parameters, d=test_doses)    
    }

    # Question- this is not created from sample dataset
    # Is it even possible to create bmd density plot for each model?
    
    temp <- fit$bmd_dist[,1]
    temp <- temp[!is.infinite(temp)]
    
    temp_density<-data.frame(matrix(0,length(temp),3))
    temp_density[,2]=fit$model
    temp_density[,1]=temp
    temp_density[,3]=A$posterior_probs[i]
    
    assign(paste("t",i,sep="_"),temp_density)
    
  }
  
  t_combine<-rbind(t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9)
  
  # This part needs to be fixed
  idx <- sample(1:9, nrow(t_1),replace=TRUE,prob=A$posterior_probs)
  
  c1<-t_1$X1
  c2<-t_2$X1
  c4<-t_3$X1
  c4<-t_4$X1
  c5<-t_5$X1
  c6<-t_6$X1
  c7<-t_7$X1
  c8<-t_8$X1
  c9<-t_9$X1
  
  
  combine_samples<-data.frame(cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9))
  
  # Select MA values
  BMD_MA<-matrix(NA,nrow(t_1),1)
  for (i in 1:nrow(t_1)){
    BMD_MA[i,1]<-combine_samples[sample(nrow(combine_samples), size=1, replace=TRUE),idx[i]]
  }
  
  BMD_MA<-data.frame(BMD_MA)
  
  t_ma<-BMD_MA %>% mutate(X1=BMD_MA,X2="Model Average",X3=1)
  #BMD_CDF should be shown here - it 
  t_ma<-t_ma %>% select(X1,X2,X3)
  
  
  
  t_combine2<-rbind(t_combine,t_ma)
  
  # From samples 
  p<-ggplot()+
    stat_density_ridges(data=t_combine2, aes(x=X1, y=fct_reorder(X2,X3,.desc=T),alpha=sqrt(X3)),
                        calc_ecdf=TRUE,quantiles=c(0.025,0.975),na.rm=T,quantile_lines=T,fill="blue")+
    xlim(c(0,quantile(t_combine$X1,0.99)))+
    geom_vline(xintercept = A$bmd[1],linetype="longdash")+
    scale_fill_manual(
      name = "Probability",
      values = c("#FF0000A0", "#A0A0A0A0", "#FF0000A0"),
      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
    )+
    labs(x="Dose Level (Dotted Line : MA BMD)",y="", title="Density plots for each fitted model (Fit type: MCMC)")+theme_classic()
  
  p2<-p+
    stat_density_ridges(data=t_combine2, aes(x=X1, y=fct_reorder(X2,X3,.desc=T), fill = factor(stat(quantile))),
                        geom = "density_ridges_gradient",
                        calc_ecdf = TRUE, quantiles = c(0.025, 0.975)
    )+ scale_fill_manual(
      name = "Probability",
      values = c("#FF0000A0", "NA", "#FF0000A0"),
      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
    )+theme(legend.position="none")
  
  
  return(p2) # Return output
  

  
}



.plot.density.BMDcontinous_MA_MCMC<-function(A){
  # Construct bmd sample plots for mcmc
  X1 <- X2 <- X3 <- NULL
  class_list <- names(A)
  fit_idx    <- grep("Individual_Model",class_list)
  qprob=0.05
  
  #Dose levels
  data<-A$Individual_Model_1$data
  doses<-data[,1]
  
  
  
  t_combine<-NA
  
  
  for (i in fit_idx){
    # Loop for the model
    fit<-A[[i]]
    test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
    #probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])
    
    
    
    if (fit$model=="hill"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_hill_f, d=test_doses)
      
    }
    if (fit$model=="exp-3"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_exp_3_f, d=test_doses)
      
    }
    if (fit$model=="exp-5"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_exp_5_f, d=test_doses)
      
    }
    if (fit$model=="power"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_power_f, d=test_doses)
      
    }
    if (fit$model=="FUNL"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.cont_FUNL_f, d=test_doses)
      
    }

    
    temp <- fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp <- temp[!is.infinite(temp)]
    
    
    
    # Q <- t(Q)
    # 
    # me <- colMeans(Q)
    # lq <- apply(Q,2,quantile, probs = qprob)
    # uq <- apply(Q,2,quantile, probs = 1-qprob)
    
    Dens =  density(temp,cut=c(max(doses)))
    # what is this 0.4 means? Scale?
    
    # normalize it?-- We don't need it actually here
    # Dens$y = Dens$y/max(Dens$y) * max(probs)
    # temp = which(Dens$x < max(doses))
    # D1_y = Dens$y[temp]
    # D1_x = Dens$x[temp]
    
    
    
    # Do I need to stack up the dataset?
    
    
    temp_density<-data.frame(matrix(0,length(temp),3))
    temp_density[,2]=substr(fit$full_model,8,999)
    temp_density[,1]=temp
    temp_density[,3]=A$posterior_probs[i]
    
    # assign(paste("t",i,sep="_"),temp_density)
    
    t<-temp_density
    
    t_combine<-rbind(t_combine,t)
    
  }

  t_combine<-t_combine[-1,]
  
  
  
  # This part is needed to get MA density plots
  # Model ~xxx..
  # Number of model should be fixed 
  
  # Grep function -- index
  fit_idx

  
  
  idx <- sample(1:length(fit_idx), length(A$Individual_Model_1$mcmc_result$BMD_samples),replace=TRUE,prob=A$posterior_probs)
  
  
  

  # rbind(A[[1]]$mcmc_result$BMD_samples)  

  df<-NA

  
  # How should I initialize this?
  for (i in 1:length(fit_idx)){
    m<-A[[i]]
    c<-m$mcmc_result$BMD_samples
    df<-data.frame(cbind(df,c))
  }
  # 
  # m1<-A$Individual_Model_1
  # c1<-m1$mcmc_result$BMD_samples
  # 
  # 
  # m2<-A$Individual_Model_2
  # c2<-m2$mcmc_result$BMD_samples
  # 
  # m3<-A$Individual_Model_3
  # c3<-m3$mcmc_result$BMD_samples
  # 
  # m4<-A$Individual_Model_4
  # c4<-m4$mcmc_result$BMD_samples
  # 
  # m5<-A$Individual_Model_5
  # c5<-m5$mcmc_result$BMD_samples

  # This part should be dynamically assigned... 
  # combine_samples2<-data.frame(cbind(c1,c2,c3,c4,c5))

  df_samples<-data.frame(df[,-1])
  
  

  
  
  # Select MA values
  
  BMD_MA<-matrix(NA,length(A$Individual_Model_1$mcmc_result$BMD_samples),1)
  
  for (i in 1:length(A$Individual_Model_1$mcmc_result$BMD_samples)){
    j<-sample(nrow(df_samples), size=1, replace=TRUE)
    BMD_MA[i,1]<-df_samples[j,idx[i]]
  }

  

  BMD_MA<-data.frame(BMD_MA)
  
  t_ma<-BMD_MA %>% mutate(X1=BMD_MA,X2="Model Average",X3=1)
  #BMD_CDF should be shown here - it 
  t_ma<-t_ma %>% select(X1,X2,X3)
  
  
  
  t_combine2<-rbind(t_combine,t_ma)
  
  
  t_combine3<-t_combine2 %>% filter(as.numeric(X3)>0.05 , as.numeric(X1)!=Inf)
  
  
  # John's comment- I want to fill the color as 
  p<-ggplot()+
    stat_density_ridges(data=t_combine3, aes(x=X1, y=fct_reorder(X2,X3,.desc=T),group=X2 ,alpha=sqrt(X3), fill=cut(X3,c(0,0.99,1))), 
                        calc_ecdf=TRUE,quantiles=c(0.025,0.975),na.rm=T,quantile_lines=T,scale=0.9)+
    xlim(c(0,quantile(t_combine$X1,0.99)))+
    geom_vline(xintercept = A$bmd[1],linetype="longdash")+
    scale_fill_manual(name="X3", values=c("(0,0.99]"="darkgrey", "(0.99,1]"="red"))+
    labs(x="Dose Level (Dotted Line : MA BMD)",y="", title="Density plots for each fitted model (Fit type: MCMC)")+
    theme_classic()
  
  
  
  p2<-p+theme(legend.position="none", axis.text.y=element_text(size=12))
  
  
  return(p2) # Return output
  
  
}


