This file contains the functions R to reproduce the results of the paper 
"Extreme sea level estimation combining systematic observed skew surges and historical record sea levels"

- model.R
  --> This code needs to be computed first to build the various stan models. 

- buildRData.R
  --> This code allows to build a RData file containing all the informations for each site and to implement
	each method.

- dataForStanModels.R
  --> This files contains the following functions : 
	- tide_empirical() --> computes the empirical distribution of astronomical high tides
	- data_gpd() --> prepares the info for the stan model corresponding to systematic skew surges
	- data_gpd_historicalSkewSurges() --> prepares the info for the stan model corresponding to systematic 
						skew surges and historical skew surges
	- data_gpd_historicalSeaLevels() --> prepares the info for the stan model corresponding to systematic 
						skew surges and historical sea levels

- methods.R
  --> This files contains the following functions : 
	- method_1(), method_2(), method_3(), method_4(), method_5(), method_6() 
	           --> simulate a sample for Monte Carlo simulations for each method and returns the 
			right stan model and the required data in the good format
	- values() --> returns the maximum likelihood estimates (parameters and 100-year quantile), 
			the 90% credibility interval of the 100-year quantile and the exceedance probability 
			of the real 100-year quantile for a Monte Carlo sample 

- validationFunctions.R
  --> This file contains the function validation() which allows to simulate 1000 Monte Carlo samples in the aim
	to to build the rank histogram and the boxplots of the ML 100-year quantile and of the ML parameters.

- validation.R
  --> This code calls the function validation() for each site and each method.

- adjustments.R
  --> This file computes the credibility intervals and the ML posterior quantiles for each site and each method.

- DIY.R
  --> This code allows to reproduce the results of the paper "Extreme sea level estimation combining systematic 
	observed skew surges and historical record sea levels"
