# --------------------------------------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file allows to build RData files containig all the information to implement each method.

# --------------------------------------------------------------------------------------------------------------------------------------

# path
setwd("")

# Choose a site ("Brest", "Dunkerque", "LaRochelle" or "SaintNazaire")
site = ""

dir_data = "./data/"

# loading the function tide_empirical()
source("./functions/dataForStanModels.R")

# --------------------------------------------------------------------------------------------------------------------------------------

# Selected parameters for the Monte Carlo simulations, 
# systematic duration, POT threshold
# historical duration, threshold for historical sea levels
x = read.table(paste(dir_data,site,"_info.txt",sep=""),sep="\t")
sigma = x$sigma ; xi = x$xi ; lambda = x$lambda
w_S = x$w_S ; u = x$u
w_H = x$w_H ; eta_H = x$eta_H ; eta_H_b = x$eta_H_b
x100 = x$x100

# Observed systematic skew surges (lower or larger than u)
x = read.table(paste(dir_data,site,"_allSkewSurge.txt",sep=""),sep="\t")
skewSurge = x$skewSurge
# Ordinary skew surges (lower than u)
surge_ = skewSurge[skewSurge<u]
p = length(surge_)/length(skewSurge)

# Predicted high tide values over a saros cycle (18,6 years)
x = read.table(paste(dir_data,site,"_maxPrediction.txt",sep=""),sep="\t")
maxPrediction = x$maxPrediction

# Empirical distribution of astronomical high tides
tide_d = tide_empirical(maxPrediction)

# Systematic skew surges larger than u
x = read.table(paste(dir_data,site,"_surgePOT.txt",sep=""),sep="\t")
skewSurgePOT = x$skewSurgePOT

# Historical information
x = read.table(paste(dir_data,site,"_HI.txt",sep=""),sep="\t")
dateH = x$date ; zH = x$maxSeaLevel ; xH = x$skewSurge

if (site=="Dunkerque" || site=="LaRochelle") {
  save(skewSurge, dateH, zH, xH, sigma, xi, lambda, w_S, u, w_H, eta_H, eta_H_b, x100, maxPrediction, skewSurgePOT, tide_d, surge_, p,
       file=paste(dir_data,site,"_data.RData",sep=""))
} else {
  # Skew surges larger than u 
  x = read.table(paste(dir_data,site,"_surgePOT_complete.txt",sep=""))
  skewSurgePOT_complete = x$surge_POT
  save(skewSurge, dateH, zH, xH, sigma, xi, lambda, w_S, u, w_H, eta_H, eta_H_b, x100, maxPrediction, skewSurgePOT, skewSurgePOT_complete, tide_d, surge_, p,
       file=paste(dir_data,site,"_data.RData",sep=""))
}
