This folder contains the data of the following sites : Brest, Dunkerque, La Rochelle and Saint Nazaire.

- site_allSkewSurge.txt
  --> all skew surges (lower and larger than u) during the systematic period

- site_HI.txt
  --> historical information (historical record sea levels and their corresponding skew surges)

- site_info.txt
  --> selected parameters (sigma, xi, lamdba) for the Monte Carlo simulations, systematic and 
	historical durations (w_S, w_H), systematic and historical sampling thresholds (u, eta_H)

- site_maxPrediction.txt 
  --> prediction of the maximum tidal level for a period of 18.6 years (a saros)

- site_surge_POT.txt
  --> POT sample of the systematic skew surges

- site_data.RData
  --> includes the data of the 5 previous files, built with the function buildRData.R in the the folder "functions"