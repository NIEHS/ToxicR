# Continuous functions are defined here

# FUNL
.cont_FUNL_f <- function(parms,doses){
  b <- parms[1] + parms[2]*exp((doses-parms[5])^2*(-parms[6]))*(1/(1+exp(-(doses-parms[3])/parms[4])))
  return(b)
}

#dichotomous hill
.cont_hill_f <- function(parms,d){
  g  <- parms[1] 
  nu <- parms[2]  
  
  k  <- parms[3];
  n  <- parms[4]; 
  rval <- g + nu*d^n/(k^n+d^n)
  return (rval)
}

#dichotomous log-logistic
.cont_exp_5_f <- function(parms,d){
  g <- parms[1]
  b <- parms[2];
  c <- parms[3];
  e <- parms[4]; 
  rval <- g*(exp(c)-(exp(c)-1.0)*(exp(-(b*d)^e)))
  return (rval)
}

#dichotomous log-probit
.cont_exp_3_f <-function(parms,d){
  g <- parms[1]
  b <- parms[2]
  e <- parms[4] 
  rval <- g*exp((b*d)^e)
  return (rval)
}

.cont_power_f <-function(parms,d){
  g <- parms[1]; 
  b <- parms[2];
  a <- parms[3]; 
  rval <- g + b*d^a
  return (rval)
}
