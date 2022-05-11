# #############################################################
# # Function: clean_continuous_analysis
# # By default a continuous analysis constraints the result
# # to be proportional to background mean.  This is 'undone'
# # in this function so all values are on the origional scale. 
# # info: the fit information for a continuous model
# # log_normal: TRUE  if the fit was a lognormal
# # dmodel:   Fit information type... to be fixed. 
# #############################################################
# clean_continuous_analysis <-function(info,log_normal,dmodel){
#   
#   c("hill","exp-3","exp-5","power")
#   if (dmodel == 1) #hill
#   {
#     if (log_normal){
#       DELTA <- diag(c(info$SCALE,info$SCALE,info$dose_scale,1,1));
#       info$EST <- DELTA%*%info$EST; 
#       info$COV <- DELTA%*%info$COV%*%t(DELTA)
#       if (ncol(info$Y)==1){
#           info$F_INFO[1] =  info$F_INFO[1] - sum(log(1/info$Y)) + sum(log(info$SCALE/info$Y)) 
#       }else{
#           
#       }
#     }else{  
#       #rescale the estimates for normality
#       if (nrow(info$EST)==5) {#constant variance
#         DELTA <- diag(c(info$SCALE,info$SCALE,info$dose_scale,1,1));
#         info$EST <- DELTA%*%info$EST; tmp = info$EST[5] + 2*log(info$SCALE) 
#         info$COV <- DELTA%*%info$COV%*%t(DELTA)
#         #fix the constant portion of the likelihood
#         if (ncol(info$Y)==1){ #individual data
#           info$F_INFO[1] =  info$F_INFO[1]+log(exp(0.5*tmp*nrow(info$Y)))-log(exp(0.5*info$EST[5]*nrow(info$Y))) 
#           info$EST[5]    =  tmp
#         }else{ #summarized data
#           
#         }
#       }else{ #non-constant variance
#         DELTA <- diag(c(info$SCALE,info$SCALE,info$dose_scale,1,1,1));
#         o_est <- info$EST
#         info$EST <- DELTA%*%info$EST;  
#         info$COV <- DELTA%*%info$COV%*%t(DELTA)
#         #fix the constant portion of the likelihood
#         if (ncol(info$Y)==1){ #individual data
#           t_est = info$EST
#           mean_t  = t_est[1] + t_est[2]*(info$X^t_est[4])/(t_est[3]^t_est[4]+(info$X[,1])^t_est[4])
#           mean_o  = o_est[1] + o_est[2]*(info$X/max(info$X[,1]))^o_est[4]/(o_est[3]^o_est[4]+(info$X/max(info$X))^o_est[4])
#           var_t   = exp(info$EST[6]+(2-info$EST[5])*log(info$SCALE))*mean_t^info$EST[5]
#           var_o   = exp(info$EST[6])*mean_o^info$EST[5]
#           info$F_INFO[1] =  info$F_INFO[1]+sum(0.5*log(var_t))-sum(0.5*log(var_o))
#           info$EST[6]    =  info$EST[6]+(2-info$EST[5])*log(info$SCALE)
#         }else{ #summarized data
#           
#         }
#       }
#     }
#     
#   }   
#   if (dmodel == 2) #exp-3
#   { #lognormality is the easiest
#     if (log_normal){
#       DELTA <- diag(c(info$SCALE,(1/info$dose_scale),1,1));
#       info$EST <- DELTA%*%info$EST; 
#       info$COV <- DELTA%*%info$COV%*%t(DELTA)
#       if (ncol(info$Y)==1){ #individual data
#         info$F_INFO[1] =  info$F_INFO[1] - sum(log(1/info$Y)) + sum(log(info$SCALE/info$Y))    
#       }else{
#         
#       }  
#     }else{  
#       #rescale the estimates for normality
#         if (nrow(info$EST)==4) {#constant variance
#           DELTA <- diag(c(info$SCALE,(1/info$dose_scale),1,1));
#           info$EST <- DELTA%*%info$EST; tmp = info$EST[4] + 2*log(info$SCALE) 
#           info$COV <- DELTA%*%info$COV%*%t(DELTA)
#           #fix the constant portion of the likelihood
#           if (ncol(info$Y)==1){ #individual data
#               info$F_INFO[1] =  info$F_INFO[1]+log(exp(0.5*tmp*nrow(info$Y)))-log(exp(0.5*info$EST[4]*nrow(info$Y))) 
#               info$EST[4]    =  tmp
#           }else{ #summarized data
#             
#           }
#         }else{ #non-constant variance
#           DELTA <- diag(c(info$SCALE,(1/info$dose_scale),1,1,1));
#           o_est <- info$EST
#           info$EST <- DELTA%*%info$EST;  
#           info$COV <- DELTA%*%info$COV%*%t(DELTA)
#           #fix the constant portion of the likelihood
#           if (ncol(info$Y)==1){ #individual data
#             t_est = info$EST
#             mean_t  = t_est[1] * exp((t_est[2]*info$Y)^t_est[3])
#             mean_o  = o_est[1] * exp((o_est[2]*info$Y/max(info$Y))^o_est[3])
#             var_t   = exp(info$EST[5]+(2-info$EST[4])*log(info$SCALE))*mean_t^info$EST[4]
#             var_o   = exp(info$EST[5])*mean_o^info$EST[4]
#             info$F_INFO[1] =  info$F_INFO[1]+sum(0.5*log(var_t))-sum(0.5*log(var_o))
#             info$EST[5]    =  info$EST[5]+(2-info$EST[4])*log(info$SCALE)
#           }else{ #summarized data
#             
#           }
#         }
#     }
#     
#   }  
#   if (dmodel == 3) #exp-5
#   {
#     if (log_normal){
#       DELTA <- diag(c(info$SCALE,(1/info$dose_scale),1,1,1));
#       info$EST <- DELTA%*%info$EST; 
#       info$COV <- DELTA%*%info$COV%*%t(DELTA)
#       if (ncol(info$Y)==1){ #individual data
#          info$F_INFO[1] =  info$F_INFO[1] - sum(log(1/info$Y)) + sum(log(info$SCALE/info$Y))
#       }else{
#         
#       }
#     }  
#     else{  
#       #rescale the estimates for normality
#       if (nrow(EST)==5) {#constant variance
#         DELTA <- diag(c(info$SCALE,(1/info$dose_scale),1,1,1));
#         info$EST <- DELTA%*%info$EST; tmp = info$EST[5] + 2*log(info$SCALE) 
#         info$COV <- DELTA%*%info$COV%*%t(DELTA)
#         #fix the constant portion of the likelihood
#         if (ncol(info$Y)==1){ #individual data
#           info$F_INFO[1] =  info$F_INFO[1]+log(exp(0.5*tmp*nrow(info$Y)))-log(exp(0.5*info$EST[5]*nrow(info$Y))) 
#           info$EST[5]    =  tmp
#         }else{ #summarized data
#           
#         }
#       }else{ #non-constant variance
#         DELTA <- diag(c(info$SCALE,(1/info$dose_scale),1,1,1,1));
#         o_est <- info$EST
#         info$EST <- DELTA%*%info$EST;  
#         info$COV <- DELTA%*%info$COV%*%t(DELTA)
#         #fix the constant portion of the likelihood
#         if (ncol(info$Y)==1){ #individual data
#           t_est = info$EST
#           mean_t  = t_est[1] *(exp(t_est[3])-(exp(t_est[3])-1)*exp(-(t_est[2]*info$Y)^t_est[4]))
#           mean_o  = o_est[1] *(exp(o_est[3])-(exp(o_est[3])-1)*exp(-(o_est[2]*info$Y/max(info$Y))^o_est[4]))
#           var_t   = exp(info$EST[6]+(2-info$EST[5])*log(info$SCALE))*mean_t^info$EST[5]
#           var_o   = exp(info$EST[6])*mean_o^info$EST[5]
#           info$F_INFO[1] =  info$F_INFO[1]+sum(0.5*log(var_t))-sum(0.5*log(var_o))
#           info$EST[6]    =  info$EST[6]+(2-info$EST[5])*log(info$SCALE)
#         }else{ #summarized data
#           
#         }
#       }
#     }
#     
#   }  
#   if (dmodel == 4) #power
#   {
#     if (log_normal){
#       DELTA <- diag(c(info$SCALE,info$SCALE*(1/info$dose_scale)^(info$EST[3]),1,1));
#       info$EST <- DELTA%*%info$EST; 
#       info$COV <- DELTA%*%info$COV%*%t(DELTA)
#       if (ncol(info$Y)==1){ #individual data
#           info$F_INFO[1] =  info$F_INFO[1] - sum(log(1/info$Y)) + sum(log(info$SCALE/info$Y))   
#       }else{
#         
#       }
#     }else{  
#       #rescale the estimates for normality
#       if (nrow(info$EST)==4) {#constant variance
#         DELTA <- diag(c(info$SCALE,info$SCALE*(1/info$dose_scale)^(info$EST[3]),1,1));
#         info$EST <- DELTA%*%info$EST; tmp = info$EST[4] + 2*log(info$SCALE) 
#         info$COV <- DELTA%*%info$COV%*%t(DELTA)
#         #fix the constant portion of the likelihood
#         if (ncol(info$data)){ #individual data
#             info$F_INFO[1] =  info$F_INFO[1]+log(exp(0.5*tmp*nrow(info$Y)))-log(exp(0.5*info$EST[4]*nrow(info$Y))) 
#             info$EST[4]    =  tmp
#         }else{ #summarized data
#           
#         }
#       }else{ #non-constant variance
#         DELTA <- diag(c(info$SCALE,info$SCALE*(1/info$dose_scale)^(info$EST[3]),1,1,1));
#         o_est <- info$EST
#         info$EST <- DELTA%*%info$EST;  
#         info$COV <- DELTA%*%info$COV%*%t(DELTA)
#         #fix the constant portion of the likelihood
#         if (ncol(info$Y)==1){ #individual data
#             t_est = info$EST
#             mean_t  = t_est[1] + t_est[2]*info$Y^t_est[3]
#             mean_o  = o_est[1] + o_est[2]*(info$Y/max(info$Y))^o_est[3]
#             var_t   = exp(info$EST[5]+(2-info$EST[4])*log(info$SCALE))*mean_t^info$EST[4]
#             var_o   = exp(info$EST[5])*mean_o^info$EST[4]
#             info$F_INFO[1] =  info$F_INFO[1]+sum(0.5*log(var_t))-sum(0.5*log(var_o))
#             info$EST[5]    =  info$EST[5]+(2-info$EST[4])*log(info$SCALE)
#         }else{ #summarized data
#           
#         }
#       }
#     }
#   }
#   
#   return(info)
# }
# 
# #############################################################
# # Function: clean_continuous_analysis
# # By default a continuous analysis constraints the result
# # to be proportional to background mean.  This is 'undone'
# # in this function so all values are on the origional scale. 
# # info: the fit information for a continuous model
# # log_normal: TRUE  if the fit was a lognormal
# # dmodel:   Fit information type... to be fixed. 
# #############################################################
# clean_parameters <-function(model,A,B,SCALE,dose_scale,deg,log_normal){
#   EST <- A 
#   COV <- B
# 
#   dmodel = which(model==c("hill","exp-3","exp-5","power","poly"))
#   
#   if (dmodel == 1) #hill
#   {   #rescale the estimates for normality
#       if (nrow(EST)==5) {#constant variance
#         DELTA <- diag(c(SCALE,SCALE,dose_scale,1,1));
#         EST <- DELTA%*%EST;  
#         if (!log_normal){
#           EST[5] = EST[5] + 2*log(SCALE)
#         }
#         COV <- DELTA%*%COV%*%t(DELTA)
#       }else{ #non-constant variance
#         DELTA <- diag(c(SCALE,SCALE,dose_scale,1,1,1));
#         EST <- DELTA%*%EST;  
#         COV <- DELTA%*%COV%*%t(DELTA)
#       }
#   }
#   if (dmodel == 2) #exp-3
#   { #lognormality is the easiest
#     if (nrow(EST)==4) {#constant variance
#         DELTA <- diag(c(SCALE,(1/dose_scale),1,1));
#         EST <- DELTA%*%EST; 
#         if (!log_normal){
#           EST[4] = EST[4] + 2*log(SCALE)
#         }
#         COV <- DELTA%*%COV%*%t(DELTA)
#     }else{ #non-constant variance
#         DELTA <- diag(c(SCALE,(1/dose_scale),1,1,1));
#         EST <- DELTA%*%EST;  
#         COV <- DELTA%*%COV%*%t(DELTA)
#     }
#   }
#   if (dmodel == 3) #exp-5
#   {  
#      if (nrow(EST)==5) {#constant variance
#         DELTA <- diag(c(SCALE,(1/dose_scale),1,1,1));
#         EST <- DELTA%*%EST; 
#         if (!log_normal){
#           EST[5] = EST[5] + 2*log(SCALE)
#         }
#         COV <- DELTA%*%COV%*%t(DELTA)
#       }else{ #non-constant variance
#         DELTA <- diag(c(SCALE,(1/dose_scale),1,1,1,1));
#         EST <- DELTA%*%EST;  
#         COV <- DELTA%*%COV%*%t(DELTA)
#       }
#   }
#   if (dmodel == 4) #power
#   {
#     if (nrow(EST) == 4) {#constant variance
#       DELTA <- diag(c(SCALE,SCALE*(1/dose_scale)^(EST[3]),1,1));
#       EST <- DELTA%*%EST;  
#       if (!log_normal){
#          EST[4] = EST[4] + 2*log(SCALE)
#       }
#       DELTA[3,4] = log(1/dose_scale)*SCALE*(1/dose_scale)^EST[3] #delta method
#       COV <- DELTA%*%COV%*%t(DELTA)
#     } else{ #non-constant variance
#       DELTA <- diag(c(SCALE,SCALE*(1/dose_scale)^(EST[3]),1,1,1));
#       EST <- DELTA%*%EST;  
#       DELTA[3,4] = log(1/dose_scale)*SCALE*(1/dose_scale)^EST[3] #delta method
#       COV <- DELTA%*%COV%*%t(DELTA)
#     }
#   }
#   if (dmodel == 5) #polynomial
#   {
#     if (nrow(EST)==deg+2) {#constant variance
#       DELTA <- diag(c(SCALE,SCALE*(1/dose_scale)^(1:deg),1));
#       EST <- DELTA%*%EST;  
#       if (!log_normal){
#         EST[4] = EST[4] + 2*log(SCALE)
#       }
#       COV <- DELTA%*%COV%*%t(DELTA)
#     } else{ #non-constant variance
#       DELTA <- diag(c(SCALE,SCALE*(1/dose_scale)^(1:deg),1,1));
#       EST <- DELTA%*%EST;  
#       COV <- DELTA%*%COV%*%t(DELTA)
#     }
#   }
#   rval <- list(EST,COV)
#   names(rval) <- c("EST","COV")
#   return(rval)
# }
