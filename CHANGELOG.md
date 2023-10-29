# Changes
## Version 23.10.1.1.1
## The following cumulative fixes are in version 23.10.1.1.1
-  The function summary() would sometimes return a different BMD/BMDL than ma_continuous_fit(),
   this is now fixed. 
-  When the profile likelihood is used and the profile did not properly bound the BMDL, the function
   could return a negative BMDL.  This behavior has been fixed, and it now reports suspect model fits
   of this form.
-  For edge datasets, the MA may return an invalid result due to difficulty in computing the model 
   weights using the Laplace approximation.  In these cases, the model weight is set to NA and removed 
   from the average.
## New Features
 - By default, ToxicR now uses 2 threads to comply with CRAN's policy. All functions allow the user to 
   specify more threads if available. This new default will slow down MA computations, but it is
   necessary to put ToxicR back on CRAN. 
## Version 23.4.1.1.0
## The following cumulative fixes are in version 23.4.1.1.0
 - Fixed log-probit fit, which was giving log-logistic plots. 
 - Fixe MA summary.  The Proper BMD (BMDL, BMDU) are now reported
## New Features
 - EFSA continuous models added with model averaging.  These are now the default
   set of models. 
 - Increased reliability of initialization.  
## Version 23.1.1.0.8
### The following cumulative fixes are in version 23.1.1.0.8
	- Removed depricated C sprintf files. 
	- Fixed ASAN memory error in initialization routine under CLANG 15.0.6
	- Removed tidyverse as a dependency. 
## Version 22.8.1.0.4
### The following fixes are in version 1.0.4
	- Fixed valgrind initialization issue with polyK test. 
	- There is no added functionality in this release.  Keeping the 22.8 moniker. 

## Version 22.8.1.0.3
### The following fixes are in version 1.0.3
	- Fixed an NLOPT warning with CLANG 15 to keep on CRAN. 
	- There is no added functionality in this release.  Keeping the 22.8 moniker. 

## Version 22.8.1.0.2

### The following fixes are in version 1.0.2:
 	- Function 'single_continuous_fit' and 'ma_continuous_fit' changed error when defining default priors
		 for 'distribution=normal-ncv' when data are negative. Originally the variance was described as mean(Y)/var(Y); 
 		however, for negative means, this causes NA error. It is now defined as abs(mean(Y))/var(Y). 
 	- Log-normal distribution fits were incorrect when summarized data was used. The correct transformation of
		summarized data is now performed. The formula for standard deviation was typed in as sqrt(log((sd)^2/mu + 1)) it is now sqrt(log((sd/mu)^2+1)). 
	- Changed default priors for dichotomous fits to be consistant with Wheeler et al. (2020). 

### The following changes to fitting were made: 
	- Changed MLE Polynomial fit behavior.  Now the terms up to the quadratic are constrained to be in the direction 
	  of the response.  After this, i.e., degree >= 3, the parameters are unconstrained. 
	- Added summary and print methods for mcmc model averaging. 

### Known Problems not yet fixed
	- GoF for MA individual models not given. 
	- GoF for dichotomous models with (0/1) data fails. 	

## Version 22.5.1.0.1

### The following fixes are in version 22.5.1.0.1:

	- Function `single_continuous_fit' fixed prior issue with Log-Normal data, when sufficient statistics are given.
	- Log-Normal deviance for Exponential 3/5 was producing incorrect values. Now reporting correct values. 
	- Function `single_dichotomous_fit' did not return bmd_dist as an element of the return object when fit_type = 'mcmc'.
	- Dichotomous MA individual models were mislabled.  They now are consistant with Continuous model averaging using the 
	  `Individual_model' naming. 
	- Inf, NA, NaN rows are now removed before fitting. 

### The following changes to fitting were made: 

	- Changed the profile likelihood stopping criteria for profile likelihood equality constrained optimization for continuous models to be 5e-5 for the absolute change in the parameter value, i.e., |\beta_n - \beta_{n+1}| < 5e-5  is the stopping criteria. 
	- When OpenMP is enabled, the fitting of single continuous models and deviance models is done with multiple threads. 
	- MCMC Metropolis-Hastings algorithm proposal changed.  Proposals are now 125% MAP variance. This change improves mixing and increases effective
	sample size by around 10 to 20%.

### The following additional functionality was added:  

	- Added summary/print functions for Single Continuous and Single Dichotomous Models.
	- Added summary/print function for Model Averaging.
	- ntp tests (Williams/Shirley/PolyK/Dunn/Dunnett) e.g. ?ntp_williams
	- Data set from Wheeler, Bailer, Piegorsh (2019) added. 

