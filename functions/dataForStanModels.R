# --------------------------------------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file contains the following functions : tide_empirical, data_gpd, data_gpd_historicalSkewSurges, data_gpd_historicalSeaLevels

# --------------------------------------------------------------------------------------------------------------------------------------

# loading of stan models
load("./stanModel/model.RData")

# --------------------------------------------------------------------------------------------------------------------------------------

tide_empirical = function(tide) {
  
  # This function returns the empirical distribution of astronomical high tides.
  
  # Input:
  #   - tide: vector of astronomical high tides predicted on a 18,6-years period
  # Output:
  #   - tide_d: matrix describing the empirical distribution of astronomical high tides by intervals of length 0.01m
  
  # definition of the intervals such as the widths of the intervals is equal to 0.01m
  tide_intervals = seq(floor(min(tide)*100)/100, ceiling(max(tide)*100)/100, by=0.01)
  n_T = length(tide_intervals) # Number of intervals
  tide_prob = seq(0,0,length.out = n_T)
  for (i in 2:n_T) {
    tide_prob[i] = length(tide[tide<tide_intervals[i]]) - sum(tide_prob[1:i-1]) # number of values (from tide) in the interval i 
  }
  tide_intervals = tide_intervals - 0.005 # median value of each interval
  tide_d = matrix(NA,nrow=n_T,ncol=2)
  tide_d[,1] = tide_intervals
  tide_d[,2] = tide_prob/length(tide) # tide_prob is divided by the total number of values (from tide)
  
  return(tide_d)
}

data_gpd = function (xS, u, w_S) {
  
  # This function returns i) the stan model corresponding to systematic skew surges and 
  #                       ii) the data in the right format for this model. 
  
  # Input:
  #   - xS: vector of systematic skew surges
  #   - u: POT threshold
  #   - w_S: systematic duration
  # Output:
  #   - data: list of data for the stan model of systematic skew surges
  #   - model: stan model of systematic skew surges
  
  data = list(u=u,n=length(xS),xS=xS,w_S=w_S)
  model = gpd_model
  
  return (list(data,model))
}

data_gpd_historicalSkewSurges = function (xS, u, w_S, u_H, w_H, xH=NULL) {
  
  # This function returns i) the stan model corresponding to systematic and historical skew surges and 
  #                       ii) the data in the right format for this model. 
  
  # Input:
  #   - xS: vector of systematic skew surges
  #   - u: POT threshold
  #   - w_S: systematic duration
  #   - u_H: threshold of the historical skew surges
  #   - w_H: historical duration
  # Output:
  #   - data: list of data for the stan model of systematic and historical skew surges
  #   - model: stan model of systematic and historical skew surges
  
  # Warning : different models according to the number of historical skew surges
  if (length(xH)>1) {
    data = list(u=u,n=length(xS),xS=xS,w_S=w_S,u_H=u_H,w_H=w_H,h_x=length(xH),xH=xH)
    model = gpd_hSS_hx_model
  } 
  else if (length(xH)==1) {
    data = list(u=u,n=length(xS),xS=xS,w_S=w_S,u_H=u_H,w_H=w_H,xH=xH)
    model = gpd_hSS_1_model
  }
  else {
    data = list(u=u,n=length(xS),xS=xS,w_S=w_S,u_H=u_H,w_H=w_H)
    model = gpd_hSS_0_model
  }
  
  return (list(data,model))
}

data_gpd_historicalSeaLevels = function(xS, u, w_S, surge_, tide_d, p, w_H, zH=NULL, eta_H) {

  # This function returns i) the stan model corresponding to systematic skew surges and historical sea levels and 
  #                       ii) the data in the right format for this model. 
  
  # Input:
  #   - xS: vector of systematic skew surges
  #   - u: POT threshold
  #   - w_S: systematic duration
  #   - surge_: vector of observed skew surges lower than u
  #   - tide_d: probability and value of tide intervals (output of the function tide_empirical)
  #   - p: probability that a skew surge is lower than u
  #   - w_H: historical duration
  #   - z_H: historical sea levels
  #   - eta_H: threshold for the historical sea levels
  # Output:
  #   - data: list of data for the stan model of systematic and historical skew surges
  #   - model: stan model of systematic and historical skew surges
  
  N = 706*w_H # number of historical sea levels (lower and larger than eta_H) during the historical period
  nn = dim(tide_d)[1]
  
  eta_H_less_tide_d = c() # historical threshold (eta_H) - tide
  P_eta_H = 0 # probability of non exceedance of eta_H when the skew surge is lower than u
  for (i in 1:nn) {
    eta_H_less_tide_d[i] = eta_H - tide_d[i,1]
    P_eta_H = P_eta_H + tide_d[i,2]*length(surge_[surge_<=eta_H_less_tide_d[i]])/length(surge_)
  }
  
  h_z = length(zH) # number of historical sea levels larger than eta_H during the historical period 
  if (h_z>=1) {
    
    xH_low = matrix(NA,nrow=nn,ncol=h_z) # historical sea levels (larger than eta_H) - tide
    xH_up = matrix(NA,nrow=nn,ncol=h_z)  # [historical sea levels (larger than eta_H)]*1.01 - tide
    P_xH_low = seq(0,0,length.out=h_z)   # probability of non exceedance of xH_low
    P_xH_up = seq(0,0,length.out=h_z)    # probability of non exceedance of xH_up
    
    for (i in 1:nn) {
      for (j in 1:h_z) {
        xH_low[i,j] = zH[j] - tide_d[i,1]
        xH_up[i,j]  = zH[j]*1.01 - tide_d[i,1]
        P_xH_low[j] = P_xH_low[j] + tide_d[i,2]*length(surge_[surge_<=xH_low[i,j]])/length(surge_)
        P_xH_up[j]  = P_xH_up[j]  + tide_d[i,2]*length(surge_[surge_<=xH_up[i,j]])/length(surge_)
      }
    }
    
    # Warning : different models according to the number of historical sea levels
    if (h_z>1) {
      data = list(u=u, n=length(xS), xS=xS, w_S=w_S,
                  nn=nn, tide_d=tide_d, p=p,
                  eta_H_less_tide_d=eta_H_less_tide_d, P_eta_H=P_eta_H, N=(w_H*706), h_z=h_z,
                  xH_low=xH_low, xH_up=xH_up, P_xH_low=P_xH_low, P_xH_up=P_xH_up)
      model = gpd_hSL_hz_model
    } 
    else {
      data = list(u=u, n=length(xS), xS=xS, w_S=w_S,
                  nn=nn, tide_d=tide_d, p=p,
                  eta_H_less_tide_d=eta_H_less_tide_d, P_eta_H=P_eta_H, N=(w_H*706),
                  xH_low=xH_low[,1], xH_up=xH_up[,1], P_xH_low=P_xH_low[1], P_xH_up=P_xH_up[1])
      model = gpd_hSL_1_model
    }
    
  } 
  else {
    data = list(u=u, n=length(xS), xS=xS, w_S=w_S,
                nn=nn, tide_d=tide_d, p=p,
                eta_H_less_tide_d=eta_H_less_tide_d, P_eta_H=P_eta_H, N=(w_H*706))
    model = gpd_hSL_0_model
  }
  
  return (list(data,model))
}
