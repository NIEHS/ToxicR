# Changes

## Version 22.04 (1.0.1)

The following bug fixes are in version 1.0.1:

	- Function `single_continuous_fit' fixed prior issue with Log-Normal data, when sufficient statistics are given.
	- Log-Normal deviance for Exponential 3/5 was producing incorrect values. Now reporting 
	correct values. 
The following changes to fitting were made: 
	- Changed the profile likelihood stopping criteria for profile likelihood equality constrained optimization to be 1e-5. 
The following additional functionality was added:  

	
	- Added summary/print functions for Single Continuous and Single Dichotomous Models.

