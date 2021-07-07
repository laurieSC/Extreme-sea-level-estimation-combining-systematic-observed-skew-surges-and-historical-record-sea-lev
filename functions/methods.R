# --------------------------------------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file contains the following functions : method_1, method_2, method_3, method_4, method_5, method_6 and values

# --------------------------------------------------------------------------------------------------------------------------------------

# package loading
library(evd)

# --------------------------------------------------------------------------------------------------------------------------------------

method_1 = function (params) {
  
  # This function simulates a sample for Monte Carlo simulations in the case of method 1 and returns the right stan model and 
  # the required data in the good format.  
  
  # Input:
  #   - params: list of following parameters 
  #     --> u: POT threshold
  #     --> lambda: Poisson parameter
  #     --> sigma: scale parameter
  #     --> xi: shape parameter
  #     --> w_S: systematic duration
  # Output:
  #   - data: data for the stan model corresponding to the method 1 
  #   - model: stan model for the implementation of method 1
  
  u = params[[1]] ; lambda = params[[2]] ; sigma = params[[3]] ; xi = params[[4]] ; w_S = params[[5]]
  
  n = rpois(1,lambda*w_S)
  xS = rgpd(n, u, sigma, xi)
  
  data_model = data_gpd(xS, u, w_S)
  data = data_model[[1]] ; model = data_model[[2]]
  
  return (list(data, model))
  
} 

method_2 = function (params) {
  
  # This function simulates a sample for Monte Carlo simulations in the case of method 2 and returns the right stan model and 
  # the required data in the good format.  
  
  # Input:
  #   - params: list of following parameters 
  #     --> u: POT threshold
  #     --> lambda: Poisson parameter
  #     --> sigma: scale parameter
  #     --> xi: shape parameter
  #     --> w_S: systematic duration
  #     --> w_H: historical duration
  # Output:
  #   - data: data for the stan model corresponding to the method 2 
  #   - model: stan model for the implementation of method 2
  
  u = params[[1]] ; lambda = params[[2]] ; sigma = params[[3]] ; xi = params[[4]] ; w_S = params[[5]] ; w_H = params[[6]]
  
  n = rpois(1,lambda*w_S)
  xS = rgpd(n, u, sigma, xi)
  h_x = rpois(1,lambda*w_H)
  xH = rgpd(h_x, u, sigma, xi)
  
  data_model = data_gpd(c(xS,xH), u, (w_S + w_H))
  data = data_model[[1]] ; model = data_model[[2]]
  
  return (list(data, model))

} 

method_3 = function (params) {
  
  # This function simulates a sample for Monte Carlo simulations in the case of method 3 and returns the right stan model and 
  # the required data in the good format.  
  
  # Input:
  #   - params: list of following parameters 
  #     --> u: POT threshold
  #     --> lambda: Poisson parameter
  #     --> sigma: scale parameter
  #     --> xi: shape parameter
  #     --> w_S: systematic duration
  #     --> w_H: historical duration
  #     --> eta_H: threshold of historical sea levels
  #     --> surge_: vector of observed skew surges lower than u
  #     --> tide: vector of astronomical high tides predicted on a 18,6-years period
  #     --> tide_d: probability and value of tide intervals (output of the function tide_empirical)
  # Output:
  #   - data: data for the stan model corresponding to the method 3
  #   - model: stan model for the implementation of method 3
  
  u = params[[1]] ; lambda = params[[2]] ; sigma = params[[3]] ; xi = params[[4]] ; w_S = params[[5]] ; w_H = params[[6]] ; eta_H = params[[7]]
  surge_ = params[[8]] ; tide = params[[9]] ; tide_d = params[[10]]
  
  n = rpois(1,lambda*w_S)
  xS = rgpd(n, u, sigma, xi)
  h_x = rpois(1,lambda*w_H)
  xH = c(sample(surge_, (706*w_H - h_x), replace=T), rgpd(h_x, u, sigma, xi))
  yH = sample(tide, 706*w_H, replace=T)
  zH = xH + yH
  ii = which(zH>=eta_H)
  zH = zH[ii]
  
  data_model = data_gpd_historicalSeaLevels(xS, u, w_S, surge_, tide_d, p, w_H, zH, eta_H)
  data = data_model[[1]] ; model = data_model[[2]]
  
  return (list(data, model))
  
} 

method_4 = function (params) {
  
  # This function simulates a sample for Monte Carlo simulations in the case of method 4 and returns the right stan model and 
  # the required data in the good format.  
  
  # Input:
  #   - params: list of following parameters 
  #     --> u: POT threshold
  #     --> lambda: Poisson parameter
  #     --> sigma: scale parameter
  #     --> xi: shape parameter
  #     --> w_S: systematic duration
  #     --> w_H: historical duration
  #     --> eta_H: threshold of historical sea levels
  #     --> surge_: vector of observed skew surges lower than u
  #     --> tide: vector of astronomical high tides predicted on a 18,6-years period
  # Output:
  #   - data: data for the stan model corresponding to the method 4
  #   - model: stan model for the implementation of method 4
  
  u = params[[1]] ; lambda = params[[2]] ; sigma = params[[3]] ; xi = params[[4]] ; w_S = params[[5]] ; w_H = params[[6]] ; eta_H = params[[7]]
  surge_ = params[[8]] ; tide = params[[9]]
  
  n = rpois(1,lambda*w_S)
  xS = rgpd(n, u, sigma, xi)
  h_x = rpois(1,lambda*w_H)
  xH = c(sample(surge_, (706*w_H - h_x), replace=T), rgpd(h_x, u, sigma, xi))
  yH = sample(tide, 706*w_H, replace=T)
  zH = xH + yH
  ii = which(zH>=eta_H)
  xH = xH[ii] ; xH = xH[xH>=u] 
  
  data_model = data_gpd_historicalSkewSurges(xS, u, w_S, min_xH, w_H, xH)
  data = data_model[[1]] ; model = data_model[[2]]
  
  return(list(data, model))
  
} 

method_5 = function (params) {
  
  # This function simulates a sample for Monte Carlo simulations in the case of method 5 and returns the right stan model and 
  # the required data in the good format.  
  
  # Input:
  #   - params: list of following parameters 
  #     --> u: POT threshold
  #     --> lambda: Poisson parameter
  #     --> sigma: scale parameter
  #     --> xi: shape parameter
  #     --> w_S: systematic duration
  #     --> w_H: historical duration
  #     --> eta_H: threshold of historical sea levels
  #     --> surge_: vector of observed skew surges lower than u
  #     --> tide: vector of astronomical high tides predicted on a 18,6-years period
  # Output:
  #   - data: data for the stan model corresponding to the method 5
  #   - model: stan model for the implementation of method 5
  
  u = params[[1]] ; lambda = params[[2]] ; sigma = params[[3]] ; xi = params[[4]] ; w_S = params[[5]] ; w_H = params[[6]] ; eta_H = params[[7]]
  surge_ = params[[8]] ; tide = params[[9]]
  
  n = rpois(1,lambda*w_S)
  xS = rgpd(n, u, sigma, xi)
  h_x = rpois(1,lambda*w_H)
  xH = c(sample(surge_, (706*w_H - h_x), replace=T), rgpd(h_x, u, sigma, xi))
  yH = sample(tide, 706*w_H, replace=T)
  zH = xH + yH
  ii = which(zH>=eta_H)
  xH = xH[ii] ; xH = xH[xH>=u] 
  
  n_w_S = length(xS)/w_S
  new_w_H = length(xH)/n_w_S
  
  data_model = data_gpd(c(xS,xH), u, (w_S + new_w_H))
  data = data_model[[1]] ; model = data_model[[2]]
  
  return(list(data, model))
  
} 

method_6 = function (params) {
  
  # This function simulates a sample for Monte Carlo simulations in the case of method 6 and returns the right stan model and 
  # the required data in the good format.  

  # Input:
  #   - params: list of following parameters 
  #     --> u: POT threshold
  #     --> lambda: Poisson parameter
  #     --> sigma: scale parameter
  #     --> xi: shape parameter
  #     --> w_S: systematic duration
  #     --> w_H: historical duration
  #     --> eta_H: threshold of historical sea levels
  #     --> surge_: vector of observed skew surges lower than u
  #     --> tide: vector of astronomical high tides predicted on a 18,6-years period
  # Output:
  #   - data: data for the stan model corresponding to the method 6
  #   - model: stan model for the implementation of method 6
  
  u = params[[1]] ; lambda = params[[2]] ; sigma = params[[3]] ; xi = params[[4]] ; w_S = params[[5]] ; w_H = params[[6]] ; eta_H = params[[7]]
  surge_ = params[[8]] ; tide = params[[9]]
  
  n = rpois(1,lambda*w_S)
  xS = rgpd(n, u, sigma, xi)
  h_x = rpois(1,lambda*w_H)
  xH = c(sample(surge_, (706*w_H - h_x), replace=T), rgpd(h_x, u, sigma, xi))
  yH = sample(tide, 706*w_H, replace=T)
  zH = xH + yH
  ii = which(zH>=eta_H)
  xH = xH[ii] ; xH = xH[xH>=u] 
  
  nb_sampled_xH = length(xH)
  if (nb_sampled_xH>=1) {
    min_xH = min(xH)
    
    fit = gpdFit(xS, u, type="mle")
    xi_fit = fit@fit$par.ests[1]
    sigma_fit = fit@fit$par.ests[2]
    p_u_H = pgpd(min_xH, u, sigma_fit, xi_fit)
    
    n_w_S = length(xS)/w_S
    p_min_xH = n_w_S*(1-p_u_H)
    new_w_H = length(xH)/p_min_xH
    
  } else {
    min_xH = u
    new_w_H = 0
    p_u_H = 0
  }
  
  while(p_u_H==1) {
  
    n = rpois(1,lambda*w_S)
    xS = rgpd(n, u, sigma, xi)
    h_x = rpois(1,lambda*w_H)
    xH = c(sample(surge_, (706*w_H - h_x), replace=T), rgpd(h_x, u, sigma, xi))
    yH = sample(tide, 706*w_H, replace=T)
    zH = xH + yH
    ii = which(zH>=eta_H)
    xH = xH[ii] ; xH = xH[xH>=u] 
    
    nb_sampled_xH = length(xH)
    if (nb_sampled_xH>=1) {
      min_xH = min(xH)
    
      fit = gpdFit(xS, u, type="mle")
      xi_fit = fit@fit$par.ests[1]
      sigma_fit = fit@fit$par.ests[2]
      p_u_H = pgpd(min_xH, u, sigma_fit, xi_fit)
      
      n_w_S = length(xS)/w_S
  
      p_min_xH = n_w_S*(1-p_u_H)
      new_w_H = length(xH)/p_min_xH
    
    } else {
      min_xH = u
      new_w_H = 0
      p_u_H = 0
    }
  }
  
  if (nb_sampled_xH>=1) {
    data_model = data_gpd_historicalSkewSurges(xS, u, w_S, min_xH, new_w_H, xH)
  } else {
    data_model = data_gpd(xS, u, w_S)
  }
  data = data_model[[1]] ; model = data_model[[2]]
  
  return(list(data, model))
  
} 

values <- function(u, sigma_, xi_, lambda_, lp_, x100) {
  
  # This function has to be called after running a stan model to have the maximum likelihood estimates (parameters and 100-year quantile), 
  # the 90% credibility interval of the 100-year quantile and the exceedance probability of the real 100-year quantile.
  
  # Input:
  #   - u: POT threshold
  #   - sigma_: posterior distribution of the scale parameter
  #   - xi_: posterior distribution of the shape parameter
  #   - lambda_: posterior distribution of the Poisson parameter
  #   - lp_: posterior distribution of the log-likelihood
  #   - x100: real 100-year quantile
  # Output:
  #   - p__: exceedance probability of the real 100-year quantile
  #   - x05_: lower bound of the 90% credibility interval of the 100-year quantile
  #   - xML_: maximum likelihood 100-year quantile
  #   - x95_: upper bound of the 90% credibility interval of the 100-year quantile
  #   - sigmaML_: maximum likelihood scale parameter
  #   - xiML_: maximum likelihood shape parameter
  #   - lambdaML_: maximum likelihood Poisson parameter
  
  # maximum likelihood estimates: parameters and 100-year quantile
  iML = which.max(lp_)
  sigmaML_ = sigma_[iML] ; xiML_ = xi_[iML] ; lambdaML_ = lambda_[iML]
  xML_ =  qgpd(1-1/(lambdaML_*100), u, sigmaML_, xiML_)
  
  # posterior distribution of the 100-year quantile
  x100_posterior = NULL
  iter_chains = length(sigma_)
  for (i in 1:iter_chains) {
    x100_posterior[i] = qgpd(1-1/(lambda_[i]*100), u, sigma_[i], xi_[i])
  }
  
  # lower and uper bounds of the 90% credibility interval of the 100-year quantile
  x05_ = quantile(x100_posterior,0.05) ; x95_ = quantile(x100_posterior,0.95)
  
  # rank of the real quantile 
  ranks = sort(c(x100_posterior,x100))
  tt = which(ranks==x100) ; gg = sample(length(tt)) ; tt = tt[gg]
  p__ = tt/length(ranks) # exceedance probability of the real quantile
  
  return(list(p__, 
              x05_, xML_, x95_,
              lambdaML_, sigmaML_, xiML_))
}
