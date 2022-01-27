#########################################
#
#
##########################################
ma_continuous_list<-function(ml,dl){

  if (length(ml) != length(dl)){
    stop("Model List Length not the same length as distribution length list.")
  }
  check_list = ml %in% c("hill","exp-3","exp-5","power","FUNL")
  if (sum(check_list) != length(ml)){
    stop('At least one model not specified as "hill","exp-3","exp-5","power",or "FUNL".')
  }
  check_dist = dl %in% c("normal","normal-ncv","lognormal")
  if (sum(check_dist) != length(dl)){
    stop('At least one distribution not in "normal","normal-ncv","lognormal".')
  }
  
  ma_list <- list()
  
  for (ii in 1:length(ml)){
    a = list(model = ml[ii], dist = dl[ii],
                prior = bayesian_prior_continuous_default(ml[ii],dl[ii]))
    class(a) <- "BMDcontinuous_bayesian_model"
    ma_list[[ii]] <- a 
  }
  
  return(ma_list)
}

print.BMDcontinuous_bayesian_model <-function(data_model){
  model_list = c("hill","exp-3","exp-5","power")
  dist_list  = c("normal","normal-ncv","lognormal")
  MODEL_NAMES = c("Hill Model","Exponential 3 Model",
                  "Exponential 5 Model","Power Model")
  DIST_NAMES = c("Normal", "Normal Non-constant Variance",
                 "Log-Normal")
  
  
  temp = data_model
  m_temp = which(temp$model  == model_list)
  d_temp = which(temp$dist   == dist_list)
  cat(sprintf("Prior for the %s\n",MODEL_NAMES[m_temp]));
  cat(sprintf("Using %s distribution.\n",DIST_NAMES[d_temp]));
  print(data_model$prior)
  
}

