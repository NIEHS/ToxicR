# Dichotomous functions are defined here
{
  .logit <- function(p)
  {
    return (log(p/(1-p)))
  }
  
  #dichotomous hill
  .dich_hill_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    n <- 1/(1+exp(-parms[2])); 
    a <- parms[3];
    b <- parms[4]; 
    rval <- g + (1-g)*n*(1/(1+exp(-a-b*log(d))))
    return (rval)
  }
  #dichotomous log-logistic
  .dich_llogist_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
    return (rval)
  }
  #dichotomous log-probit
  .dich_lprobit_f <-function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
    return (rval)
  }
  
  #dichotomous weibull
  .dich_weibull_f <-function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*(1-exp(-b*d^a))
    return (rval)
  }
  
  #dichotomous gamma
  .dich_gamma_f <-function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*pgamma(b*d,a,1)
    return (rval)
  }
  
  #dichtomous logistic
  .dich_logist_f <- function(parms,d){
    rval <- 1/(1+exp(-parms[1]-parms[2]*d))
    return (rval)
  }
  
  #dichtomous probit
  .dich_probit_f <- function(parms,d){
    rval <- pnorm(parms[1]+parms[2]*d)
    return (rval)
  }
  
  .dich_qlinear_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    return (g + (1-g)*1-exp(-a*d))
  }
  
  .dich_multistage_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    rval = d*0
    for (ii  in 2:length(parms)){
      rval = rval - parms[ii]*d^(ii-1)
    }
    return (g + (1-g)*1-exp(rval))
  }
  
  
}