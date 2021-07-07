This folder contains the various stan models :

- gpd.stan
  --> POT sample of systematic skew surges

- gpd_historicalSeaLevels_0.stan
  --> POT sample of systematic skew surges + 0 historical sea level exceeding eta_H 

- gpd_historicalSeaLevels_1.stan
  --> POT sample of systematic skew surges + 1 historical sea level exceeding eta_H 
 
- gpd_historicalSeaLevels_hz.stan
  --> POT sample of systematic skew surges + h_z historical sea levels exceeding eta_H 

- gpd_historicalSkewSurges_0.stan
  --> POT sample of systematic skew surges + 0 historical skew surge exceeding u_H 

- gpd_historicalSkewSurges_1.stan
  --> POT sample of systematic skew surges + 1 historical skew surge exceeding u_H  

- gpd_historicalSkewSurges_hx.stan
  --> POT sample of systematic skew surges + h_x historical skew surges exceeding u_H  

- model.R
  --> this R file allow to build these 7 stan models, needs to be run before the first 
	call of the stan models