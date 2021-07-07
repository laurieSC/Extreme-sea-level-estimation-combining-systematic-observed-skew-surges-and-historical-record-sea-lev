# --------------------------------------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file compute the adjustments on the several case studies

# --------------------------------------------------------------------------------------------------------------------------------------

# Choose the site of observations and the method to implement
# Note that method 2 is only possible with Brest and Saint Nazaire
site = "" # "Brest", "Dunkerque", "La Rochelle" or "Saint Nazaire"
method = ""   # "1", "2", "3", "3b", "4", "5" or "6"

# data loading
load(paste("./data/",site,"_data.RData",sep=""))

# packages loading
library(evd)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

# loading of stan models
load("./stanModel/model.RData")

# functions loading
source("./functions/dataForStanModels.R")

dir_export = paste("./output/method",method,"/",sep="")

iterations = 100000
warmup = 50000

# --------------------------------------------------------------------------------------------------------------------------------------

# model and data 
if (method=="1") {
  data_model = data_gpd(skewSurgePOT, u, w_S)
} else if (method=="2") {
  data_model = data_gpd(skewSurgePOT_complete, u, (w_S+w_H))
} else if (method=="3") {
  data_model = data_gpd_historicalSeaLevels(skewSurgePOT, u, w_S, surge_, tide_d, p, w_H, zH, eta_H)
} else if (method=="3b") {
  data_model = data_gpd_historicalSeaLevels(skewSurgePOT, u, w_S, surge_, tide_d, p, w_H, zH=NULL, eta_H_b)
} else if (method=="4") {
  data_model = data_gpd_historicalSkewSurges(skewSurgePOT, u, w_S, min(xH[xH>=u]), w_H, xH[xH>=u])
} else if (method=="5") {
  h_x = length(xH[xH>=u])
  wH_method5 = h_x/lambda
  data_model = data_gpd(c(skewSurgePOT,xH[xH>=u]), u, c(w_S+wH_method5))
} else if (method=="6") {
  h_x = length(xH[xH>=u])
  u_H = min(xH[xH>=u])
  wH_method6 = h_x/(lambda*(1 - pgpd(u_H, u, sigma, xi)))
  data_model = data_gpd_historicalSkewSurges(skewSurgePOT, u, w_S, u_H, wH_method6, xH[xH>=u])
}

# MCMC 
stanfit = stan(fit=data_model[[2]], data=data_model[[1]], chains=4, cores=4, iter=iterations, warmup=warmup)

# posterior distributions of the parameters and the log-likelihood 
mcmc = as.matrix(stanfit)
sigma_ = mcmc[,"sigma"] ; xi_ = mcmc[,"xi"] ; lambda_ = mcmc[,"lambda"] ; lp_ = mcmc[,"lp__"]

pars = data.frame(sigma_, xi_, lambda_, lp_)
colnames(pars) = c("sigma", "xi", "lambda", "lp_")
write.table(pars, file=paste(dir_export,site,"_pars_method",method,".txt",sep=""))

# adjustments
pars = read.table(paste(dir_export,site,"_pars_method",method,".txt",sep=""))
sigma_ = pars$sigma ; xi_ = pars$xi ; lambda_ = pars$lambda ; lp_ = pars$lp
# posterior distribution 
y = c(c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9),c(2:1000))
lambda_y = lambda*y[(lambda*y)>1]
x = matrix(NA,nrow=length(sigma_),ncol=length(lambda_y))
for (i in 1:length(sigma_)) {
  x[i,] = qgpd((1-1/lambda_y),u,sigma_[i],xi_[i])
}
# posterior credibility intervals
IC = apply(x, 2, FUN=quantile, probs=c(0.05,0.95))
# ML quantiles
iML = which.max(lp_) ; sigma_ML = sigma_[iML] ; xi_ML = xi_[iML]
x_ML = x[iML,]

adjustments = data.frame(y[(lambda*y)>1],IC[1,],IC[2,],x_ML)
colnames(adjustments) = c("T","IC_low","IC_up","x_ML")
write.table(adjustments,file=paste(dir_export,site,"_adjustments_method",method,".txt",sep=""))
