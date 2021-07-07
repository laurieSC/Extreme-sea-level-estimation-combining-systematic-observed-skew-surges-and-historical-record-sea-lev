# --------------------------------------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file contains the function validation

# --------------------------------------------------------------------------------------------------------------------------------------

# packages loading
library(fExtremes)
library(evd)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

# functions loading
source("./functions/dataForStanModels.R")
source("./functions/methods.R")

dir_export = "./output/"

# --------------------------------------------------------------------------------------------------------------------------------------

validation <- function(method, site, 
                       iterations=10000, warmup=iterations/2) {
  
  # This function allow the values to build the rank histogram and the boxplots of the 100-year quantile
  # and of the parameters. 
  
  # Input:
  #   - method: method to implement ("1", "2", "3", "4", "5" or "6")
  #   - w_S: systematic duration
  #   - sigma: scale parameter
  #   - xi: shape parameter
  #   - lambda: Poisson parameter
  #   - site: site of observations ("Brest", "Dunkerque", "LaRochelle", "SaintNazaire")
  #   - iterations: total number of iterations per chain
  #   - warmup: number of warmup iterations per chain
  
  # Exceedance probability of the real 100-year quantile
  excess_proba = NULL
  
  # Bounds of the 90% credibility interval and ML of the 100-year quantile
  x05 = x95 = xML = NULL
  
  # ML parameters
  lambdaML = sigmaML = xiML = NULL
  
  # According to the method:
  # function_method --> simulates samples, selects the right stan model and put in the right format the data 
  # params --> list of parameters corresponding to the selected function_method
  if (method=="1") {
    function_method = method_1
    params = list(u, lambda, sigma, xi, w_S)
  } else if (method=="2") {
    function_method = method_2
    params = list(u, lambda, sigma, xi, w_S, w_H)
  } else if (method=="3") {
    function_method = method_3
    params = list(u, lambda, sigma, xi, w_S, w_H, eta_H, surge_, maxPrediction, tide_d)
  } else if (method=="4") {
    function_method = method_4
    params = list(u, lambda, sigma, xi, w_S, w_H, eta_H, surge_, maxPrediction)
  } else if (method=="5") {
    function_method = method_5
    params = list(u, lambda, sigma, xi, w_S, w_H, eta_H, surge_, maxPrediction)
  } else if (method=="6") {
    function_method = method_6
    params = list(u, lambda, sigma, xi, w_S, w_H, eta_H, surge_, maxPrediction)
  }

  for (cpt in 1:1000) {
    
    # model and data of the simulated samples
    info = function_method(params)
    
    # MCMC 
    stanfit = stan(fit=info[[2]], data=info[[1]], chains=4, cores=4, iter=iterations, warmup=warmup)
    # posterior distributions of the parameters and the log-likelihood 
    mcmc = as.matrix(stanfit)
    sigma_ = mcmc[,"sigma"] ; xi_ = mcmc[,"xi"] ; lambda_ = mcmc[,"lambda"] ; lp_ = mcmc[,"lp__"]
    
    # Computation of the exceedance probability of the real 100-year quantile, the bounds of the 90% credibility interval 
    # and ML of the 100-year quantile and the ML parameters. 
    val = values(u, sigma_, xi_, lambda_, lp_, x100)
    excess_proba[cpt] = val[[1]]
    x05[cpt] = val[[2]] ; xML[cpt] = val[[3]] ; x95[cpt] = val[[4]]
    lambdaML[cpt] = val[[5]] ; sigmaML[cpt] = val[[6]] ; xiML[cpt] = val[[7]]
    
  }
  
  # saving the exceedance probability of the real 100-year quantile, the bounds of the 90% credibility interval 
  # and ML of the 100-year quantile and the ML parameters. 
  res = data.frame(excess_proba, x05, xML, x95, lambdaML, sigmaML, xiML)
  write.table(res,
              paste(dir_export,"method",method,"/",site,"_valid_method",method,".txt",sep=""))
  
}
