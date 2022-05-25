#Set the default clevland_plot method generic for all of the classes. 

#' Create a Cleveland plot from a model averaged model.
#'
#' @title cleveland_plot - Create a Cleveland plot from a model averaged model.
#' @param A the model averaged model to plot
#' @return Returns a \code{ggplot2} graphics object. 
#' @examples 
#' mData <- matrix(c(0, 2,50,
#'                   1, 2,50,
#'                   3, 10, 50,
#'                   16, 18,50,
#'                   32, 18,50,
#'                   33, 17,50),nrow=6,ncol=3,byrow=TRUE)
#' D <- mData[,1]
#' Y <- mData[,2]
#' N <- mData[,3]
#' 
#' model = ma_dichotomous_fit(D,Y,N)
#' cleveland_plot(model)
#' 
#' @export
cleveland_plot <- function (A){
  UseMethod("cleveland_plot")
}

# Patch note 02/03/2021 
# 1. # of model parameter input is assigned dynamically in the plot data structure
# 2. For the continous case, assign prior information only first

.cleveland_plot.BMDdichotomous_MA <- function(A){
  # Construct bmd sample plots for mcmc 
  class_list <- names(A) 
  
  # Remove "No Visible Bindings Note"
  X1 <- X2 <- X3 <-X4 <- X5 <- NULL
  # This part should be consistent
 
  fit_idx<- grep("Individual_Model",class_list)
  
  
  # Create an empty matrix to contain BMD information from each model
  bmd_ind<-matrix(0,length(fit_idx)+1,5)
  
  for (i in 1:length(fit_idx)){
    # BMD -Median
    bmd_ind[i,1]<-A[[i]]$bmd[1]
    # BMD -5%
    bmd_ind[i,2]<-A[[i]]$bmd[2]
    # BMD -95%
    bmd_ind[i,3]<-A[[i]]$bmd[3]
    # Model name 
    bmd_ind[i,4]<-A[[i]]$model
    bmd_ind[i,5]<-A$posterior_probs[i]
  }
  
  
  # Ask Matt- For the case of Laplace- bmd object is missing
  # This part should be consistent

    
  bmd_ind[length(fit_idx)+1,1]<-A$bmd[1]
  bmd_ind[length(fit_idx)+1,2]<-A$bmd[2]
  bmd_ind[length(fit_idx)+1,3]<-A$bmd[3]
  
  bmd_ind[length(fit_idx)+1,4]<-"Model Average"
  bmd_ind[length(fit_idx)+1,5]<-1
  
  bmd_ind_df<-data.frame(bmd_ind)
  bmd_ind_df$X1
  

  #Temporarily it choose from CDF case, but this should be updated
  # 
  # 
  # 
  # bmd_ind_df1 <- 
  #   bmd_ind_df %>% filter(as.numeric(X5)>0.05) 
  # 
  # bmd_ind_df2<-bmd_ind_df1 %>% mutate(X5_numeric=as.numeric(X5))
  # 
  # step1<-ggplot(bmd_ind_df2, aes(x=as.numeric(X1),y=fct_reorder(X4,X5,.desc=T)))+
  #   geom_point(aes(colour=cut(X5_numeric,c(0,0.999,1)),size=(sqrt(as.numeric(X5)))))+
  #   scale_colour_manual(name="X5_numeric", values=c("(0,0.999]"="grey", "(0.999,1]"="red"))
  # 
  # 
  # #step1
  # 
  # step2<-step1+
  #   geom_vline(xintercept=as.numeric(A$bmd[1]), linetype="dotted")+
  #   theme_minimal()+
  #   labs(x="Dose Level",y="",title="BMD Estimates by Each Model (Sorted by Posterior Probability)",size="Posterior Probability")+
  #   theme(legend.position="none")+
  #   geom_errorbar(data=bmd_ind_df, width=0.2,aes(xmin=as.numeric(X2), xmax=as.numeric(X3), y=fct_reorder(X4,X5,.desc=T)),color="blue",alpha=0.3)
  # 
  # return(step2)
  
  bmd_ind[length(fit_idx)+1,4]<-"Model Average"
  bmd_ind[length(fit_idx)+1,5]<-1
  
  bmd_ind_df<-data.frame(bmd_ind)
  
  bmd_ind_df <- 
    bmd_ind_df %>% filter(as.numeric(X5)>0.05) 
  
  #Drop NA value
  bmd_ind_df2<-data.frame(bmd_ind_df[which(!is.na(bmd_ind_df[,1])),])
  
  out<-ggplot()+
    geom_point(data=bmd_ind_df2, aes(x=as.numeric(X1), y=fct_reorder(X4,as.numeric(X5),.desc=T),size=(sqrt(as.numeric(X5)))), color="red")+
    #scale_colour_gradient(low = "gray", high = "black")+
    geom_vline(xintercept=as.numeric(bmd_ind_df2$X1[length(fit_idx)+1]), linetype="dotted")+
    theme_minimal()+
    labs(x="Dose Level",y="", title="BMD Estimates by Each Model (Sorted by Posterior Probability)",size="Posterior Probability")+
    theme(legend.position="none")+
    geom_errorbar(data=bmd_ind_df2, width=0.2,aes(xmin=as.numeric(X2), xmax=as.numeric(X3), y=fct_reorder(X4,X5,.desc=T)),color="black",alpha=0.3)
  
  return(out)
}


# Continous Case

# There should be different assumptions-Normality, Normality-NCV, Log-Normality;
# Here let's just plot Normality assumption

.cleveland_plot.BMDcontinous_MA<-function(A){
  # Construct bmd sample plots for mcmc 
  class_list <- names(A)
  # Remove "No Visible Bindings Note"
  X1 <- X2 <- X3 <-X4 <- X5 <- NULL
  # Grap function extract # of indices from the text with same pattern
  fit_idx    <- grep("Individual_Model",class_list)
  
  # Create an empty matrix to contain BMD information from each model
  bmd_ind<-matrix(0,length(fit_idx)+1,5)
  
  for (i in 1:length(fit_idx)){
    # BMD -Median
    bmd_ind[i,1]<-A[[i]]$bmd[1]
    # BMD -5%
    bmd_ind[i,2]<-A[[i]]$bmd[2]
    # BMD -95%
    bmd_ind[i,3]<-A[[i]]$bmd[3]
    # Model name 
    bmd_ind[i,4]<-substr(A[[i]]$full_model, 8,999)      
    bmd_ind[i,5]<-A$posterior_probs[i]
  }
  
  bmd_ind[length(fit_idx)+1,1]<-A$bmd[1]
  bmd_ind[length(fit_idx)+1,2]<-A$bmd[2]
  bmd_ind[length(fit_idx)+1,3]<-A$bmd[3]
  
  # Add model average case 
  bmd_ind[length(fit_idx)+1,4]<-"Model Average"
  bmd_ind[length(fit_idx)+1,5]<-1
  
  bmd_ind_df<-data.frame(bmd_ind)
  
  bmd_ind_df <- 
    bmd_ind_df %>% filter(as.numeric(X5)>0.05) 
    
  #Drop NA value
  bmd_ind_df2<-data.frame(bmd_ind_df[which(!is.na(bmd_ind_df[,1])),])
  
  out<-ggplot()+
    geom_point(data=bmd_ind_df2, aes(x=as.numeric(X1), y=fct_reorder(X4,as.numeric(X5),.desc=T),size=(sqrt(as.numeric(X5)))), color="red")+
    #scale_colour_gradient(low = "gray", high = "black")+
    geom_vline(xintercept=as.numeric(bmd_ind_df2$X1[length(fit_idx)+1]), linetype="dotted")+
    theme_minimal()+
    labs(x="Dose Level",y="", title="BMD Estimates by Each Model (Sorted by Posterior Probability)",size="Posterior Probability")+
    theme(legend.position="none")+
    geom_errorbar(data=bmd_ind_df2, width=0.2,aes(xmin=as.numeric(X2), xmax=as.numeric(X3), y=fct_reorder(X4,X5,.desc=T)),color="black",alpha=0.3)
  
  
  
  return(out)
}
