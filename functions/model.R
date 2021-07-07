# --------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file allows to build the various stan models.

# --------------------------------------------------------------------------------------------

# path
setwd("")

setwd("./stanModel/")

# loading packages
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
library(evd)

# --------------------------------------------------------------------------------------------

# characte
w_S = 50 ; u = 0.7
sigma = 0.1 ; xi = -0.2 ; lambda = 1
w_H = 100

# --------------------------------------------------------------------------------------------

model_file = "gpd.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

data_model = list(u=u,n=n,xS=xS,w_S=w_S)
gpd_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

model_file = "gpd_historicalSkewSurges_hx.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

n2 = rpois(1,w_H*lambda)
xH = rgpd(n2, u, sigma, xi) ; xH = sort(xH)[(length(xH)-1):length(xH)]
h_x = length(xH)
u_H = min(xH)

data_model = list(u=u,n=n,xS=xS,w_S=w_S,u_H=u_H,w_H=w_H,h_x=h_x,xH=xH)
gpd_hSS_hx_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

model_file = "gpd_historicalSkewSurges_1.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

n2 = rpois(1,w_H*lambda)
xH = rgpd(n2, u, sigma, xi) ; xH = sort(xH)[length(xH)]
u_H = min(xH)

data_model = list(u=u,n=n,xS=xS,w_S=w_S,u_H=u_H,w_H=w_H,xH=xH)
gpd_hSS_1_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

model_file = "gpd_historicalSkewSurges_0.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

n2 = rpois(1,w_H*lambda)
xH = rgpd(n2, u, sigma, xi) ; xH = sort(xH)[length(xH)]
u_H = min(xH) + 0.5

data_model = list(u=u,n=n,xS=xS,w_S=w_S,u_H=u_H,w_H=w_H)
gpd_hSS_0_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

model_file = "gpd_historicalSeaLevels_hz.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

nn = 201
tide_d = matrix(NA,nrow=nn,ncol=2) ; tide_d[,1] = seq(4.5,6.5,by=0.01) ; tide_d[,2] = seq(1,1,length.out=nn)/nn

p = 705/706
zS_low = sample(tide_d[,1],705*10000,replace=T) + runif(705*10000,-2,u)

n2 = rpois(1,w_H*lambda)
xH = c(runif((706*w_H - n2),-2,u), rgpd(n2, u, sigma, xi))
yH = sample(tide_d[,1], 706*w_H, replace=T)
zH = xH + yH
zH = sort(zH)[length(zH)-1:length(zH)]

eta_H = min(zH)
xH_low = matrix(NA,nrow=nn,ncol=2) ; xH_up = matrix(NA,nrow=nn,ncol=2)
eta_H_less_tide_d = c()
for (i in 1:nn) {
  eta_H_less_tide_d[i] = eta_H - tide_d[i,1];
  for (j in 1:2) {
    xH_low[i,j] = zH[j]*0.99 - tide_d[i,1]
    xH_up[i,j]  = zH[j]*1.01 - tide_d[i,1]
  }
}

P_eta_H = length(zS_low[zS_low<=eta_H])/length(zS_low)

P_xH_low = c() ; P_xH_up = c()
for (i in 1:2) {
  P_xH_low[i] = length(zS_low[zS_low<=zH[i]*0.99])/length(zS_low)
  P_xH_up[i]  = length(zS_low[zS_low<=zH[i]*1.01])/length(zS_low)
}

data_model = list(u=u, n=n, xS=xS, w_S=w_S,
                nn=nn, tide_d=tide_d, p=p,
                eta_H_less_tide_d=eta_H_less_tide_d, P_eta_H=P_eta_H, N=(w_H*706), h_z=2,
                xH_low=xH_low, xH_up=xH_up, P_xH_low=P_xH_low, P_xH_up=P_xH_up)

gpd_hSL_hz_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

model_file = "gpd_historicalSeaLevels_1.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

nn = 201
tide_d = matrix(NA,nrow=nn,ncol=2) ; tide_d[,1] = seq(4.5,6.5,by=0.01) ; tide_d[,2] = seq(1,1,length.out=nn)/nn

p = 705/706
zS_low = sample(tide_d[,1],705*10000,replace=T) + runif(705*10000,-2,u)

n2 = rpois(1,w_H*lambda)
xH = c(runif((706*w_H - n2),-2,u), rgpd(n2, u, sigma, xi))
yH = sample(tide_d[,1], 706*w_H, replace=T)
zH = xH + yH
zH = sort(zH)[length(zH)]

eta_H = min(zH)
xH_low = c() ; xH_up = c()
eta_H_less_tide_d = 0
for (i in 1:nn) {
  eta_H_less_tide_d[i] = eta_H - tide_d[i,1];
  xH_low[i] = zH*0.99 - tide_d[i,1]
  xH_up[i]  = zH*1.01 - tide_d[i,1]
}

P_eta_H = length(zS_low[zS_low<=eta_H])/length(zS_low)
P_xH_low = length(zS_low[zS_low<=zH*0.99])/length(zS_low)
P_xH_up  = length(zS_low[zS_low<=zH*1.01])/length(zS_low)

data_model = list(u=u, n=n, xS=xS, w_S=w_S,
                  nn=nn, tide_d=tide_d, p=p,
                  eta_H_less_tide_d=eta_H_less_tide_d, P_eta_H=P_eta_H, N=(w_H*706),
                  xH_low=xH_low, xH_up=xH_up, P_xH_low=P_xH_low, P_xH_up=P_xH_up)

gpd_hSL_1_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

model_file = "gpd_historicalSeaLevels_0.stan"
# --------------------------------

n = rpois(1,w_S*lambda)
xS = rgpd(n, u, sigma, xi)

nn = 201
tide_d = matrix(NA,nrow=nn,ncol=2) ; tide_d[,1] = seq(4.5,6.5,by=0.01) ; tide_d[,2] = seq(1,1,length.out=nn)/nn

p = 705/706
zS_low = sample(tide_d[,1],705*10000,replace=T) + runif(705*10000,-2,u)

n2 = rpois(1,w_H*lambda)
xH = c(runif((706*w_H - n2),-2,u), rgpd(n2, u, sigma, xi))
yH = sample(tide_d[,1], 706*w_H, replace=T)
zH = xH + yH
zH = sort(zH)[length(zH)] + 0.5

eta_H = min(zH)
eta_H_less_tide_d = 0
for (i in 1:nn) {
  eta_H_less_tide_d[i] = eta_H - tide_d[i,1]
}

P_eta_H = length(zS_low[zS_low<=eta_H])/length(zS_low)

data_model = list(u=u, n=n, xS=xS, w_S=w_S,
                  nn=nn, tide_d=tide_d, p=p,
                  eta_H_less_tide_d=eta_H_less_tide_d, P_eta_H=P_eta_H, N=(w_H*706))

gpd_hSL_0_model = stan(model_file, data=data_model, chains = 0)

# --------------------------------------------------------------------------------------------

save(gpd_model, 
     gpd_hSS_hx_model, gpd_hSS_1_model, gpd_hSS_0_model,
     gpd_hSL_hz_model, gpd_hSL_1_model, gpd_hSL_0_model, 
     file="model.RData")
