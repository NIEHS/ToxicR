# Changes

## Version 22.05 (1.0.1)

### The following fixes are in version 1.0.1:

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

