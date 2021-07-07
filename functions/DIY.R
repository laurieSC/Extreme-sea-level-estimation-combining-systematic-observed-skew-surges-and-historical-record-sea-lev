# -----------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# Reproduce the figures of the paper
# Extreme sea level estimation combining systematic observed skew surges and historical record sea levels

# -----------------------------------------------------------------------------------------------------------

# packages loading
library(evd)
library(latex2exp)

# path
setwd("")

# Choose the site ("Brest", "Dunkerque", "LaRochelle", "SaintNazaire")
site = "LaRochelle"

# LOADING 
## data 
load(paste("./data/",site,"_data.RData",sep=""))
# results of the Monte Carlo simulations
valid_method1 = read.table(paste("./output/method1/",site,"_valid_method1.txt",sep=""))
valid_method2 = read.table(paste("./output/method2/",site,"_valid_method2.txt",sep=""))
valid_method3 = read.table(paste("./output/method3/",site,"_valid_method3.txt",sep=""))
valid_method4 = read.table(paste("./output/method4/",site,"_valid_method4.txt",sep=""))
valid_method5 = read.table(paste("./output/method5/",site,"_valid_method5.txt",sep=""))
valid_method6 = read.table(paste("./output/method6/",site,"_valid_method6.txt",sep=""))

# ------------------------------------------------------------------------------------------

# Table 3 : Characteristics of the historical data sets
xH = xH[zH>=eta_H]
h_x = length(xH[xH>=u])
u_H = min(xH[xH>=u])
wH_method5 = h_x/lambda
wH_method6 = h_x/(lambda*(1 - pgpd(u_H, u, sigma, xi)))

# Table 4 : Characteristics of the generated historical series of sea levels and skew surges
min_zH_ = c() ; h_z_ = c()
min_xH_ = c() ; nb_surges_u_ = c() ; nb_surges_u_H_= c() ; h_x_ = c() ; sampling_rate_ = c()
for (cpt in 1:1000) {
  n_ = rpois(1,lambda*w_S)
  xS_ = rgpd(n_, u, sigma, xi)
  nb_surges_u_[cpt] = rpois(1,lambda*w_H)
  xH_ = c(sample(surge_, (706*w_H - nb_surges_u_[cpt]), replace=TRUE), rgpd(nb_surges_u_[cpt], u, sigma, xi))
  nb_surges_u_H_[cpt] = length(xH_[xH_>=u_H])
  yH_ = sample(maxPrediction, 706*w_H, replace=TRUE)
  zH_ = xH_ + yH_
  ii = which(zH_>=eta_H)
  zH_ = zH_[ii]
  if (length(zH_)>=1) {
    min_zH_[cpt] = min(zH_)
  } else {
    min_zH_[cpt] = eta_H
  }
  h_z_[cpt] = length(zH_)
  xH_ = xH_[ii] ; xH_ = xH_[xH_>=u] 
  if (length(xH_)>=1) {
    min_xH_[cpt] = min(xH_)
  } else {
    min_xH_[cpt] = u
  }
  
  h_x_[cpt] = length(xH_[xH_>=u_H])
}
# Generated historical sea levels 
mean(min_zH_) # Minimum generated value 
mean(h_z_)    # Average number of record values
# Generated historical skew surges
mean(min_xH_) # Minimum sampled value u_H
mean(nb_surges_u_) # Average number of skew surges > u
mean(nb_surges_u_H_) # Average number of skew surges > u_H
mean(h_x_) # Average number of sampled values
mean(h_x_/nb_surges_u_H_)*100 # Average skew surge sampling rate

# Figure 3 : Dispersion of the 100-year quantile estimated with the maximum likelihood (divided by the real value), obtained from simulations.
toto = data.frame(c(valid_method1$xML/x100,valid_method2$xML/x100,valid_method3$xML/x100,valid_method4$xML/x100,valid_method5$xML/x100,valid_method6$xML/x100),
                  c(seq.int(from=1, to=1, length.out=1000),seq.int(from=2, to=2, length.out=1000),
                    seq.int(from=3, to=3, length.out=1000),seq.int(from=4, to=4, length.out=1000),
                    seq.int(from=5, to=5, length.out=1000),seq.int(from=6, to=6, length.out=1000)))
names(toto) = c("Y","X")
plot(c(1,1),c(-10,10),
     type="l",col="lightgrey",lty=3, xlim=c(0.5,6.5),ylim=c(0.5,2),xlab="Method",ylab=TeX('$\\hat{x}^{ML}_{100}/x_{100}$'),axes=FALSE)
axis(2)
axis(1,at=c(1:6),label=c("1","2","3","4","5","6"))
lines(c(2,2),c(-10,10),col="lightgrey",lty=3)
lines(c(3,3),c(-10,10),col="lightgrey",lty=3)
lines(c(4,4),c(-10,10),col="lightgrey",lty=3)
lines(c(5,5),c(-10,10),col="lightgrey",lty=3)
lines(c(6,6),c(-10,10),col="lightgrey",lty=3)
boxplot(toto$Y~toto$X,col=c("grey","violetred1","darkolivegreen1","coral","cadetblue1","darkblue"),axes=FALSE,add=TRUE)
lines(c(-1,10),c(1,1))
box()

# Figure 4 and Table D1 : Relative bias, RSD and RRMSE of the ML estimated 100-year quantile 
relative_bias = c(mean((valid_method1$xML-x100)/x100), mean((valid_method2$xML-x100)/x100), mean((valid_method3$xML-x100)/x100), 
         mean((valid_method4$xML-x100)/x100), mean((valid_method5$xML-x100)/x100), mean((valid_method6$xML-x100)/x100))
RSD = c(sd(valid_method1$xML)/mean(valid_method1$xML), sd(valid_method2$xML)/mean(valid_method2$xML), sd(valid_method3$xML)/mean(valid_method3$xML),
        sd(valid_method4$xML)/mean(valid_method4$xML), sd(valid_method5$xML)/mean(valid_method5$xML), sd(valid_method6$xML)/mean(valid_method6$xML))
RMSE <- function(pred,Q) {
  return (sqrt(sum((pred-Q)**2/length(pred))))
}
RRMSE = c(RMSE(valid_method1$xM,x100)/mean(valid_method1$xML), RMSE(valid_method2$xM,x100)/mean(valid_method2$xML), RMSE(valid_method3$xM,x100)/mean(valid_method3$xML), 
          RMSE(valid_method4$xM,x100)/mean(valid_method4$xML), RMSE(valid_method5$xM,x100)/mean(valid_method5$xML), RMSE(valid_method6$xM,x100)/mean(valid_method6$xML))
criterion = matrix(data=NA,ncol=3,nrow=6)
criterion[,1] = relative_bias
criterion[,2] = RSD
criterion[,3] = RRMSE
round(criterion,2)

# Figure 5 (E1, E2, E3) : Dispersion of the parameters estimated with the maximum likelihood (divided by the real values), obtained from simulations with different tested methods.
## Lambda 
toto = data.frame(c(valid_method1$lambdaML/lambda,valid_method2$lambdaML/lambda,valid_method3$lambdaML/lambda,valid_method4$lambdaML/lambda,valid_method5$lambdaML/lambda,valid_method6$lambdaML/lambda),
                  c(seq.int(from=1, to=1, length.out=1000),seq.int(from=2, to=2, length.out=1000),
                    seq.int(from=3, to=3, length.out=1000),seq.int(from=4, to=4, length.out=1000),
                    seq.int(from=5, to=5, length.out=1000),seq.int(from=6, to=6, length.out=1000)))
names(toto) = c("Y","X")
plot(c(1,1),c(-10,10),
     type="l",col="lightgrey",lty=3, xlim=c(0.5,6.5),ylim=c(0.5,2),xlab="Method",ylab=TeX('$\\hat{\\lambda}^{ML}/\\lambda$'),axes=FALSE)
axis(2)
axis(1,at=c(1:6),label=c("1","2","3","4","5","6"))
lines(c(2,2),c(-10,10),col="lightgrey",lty=3)
lines(c(3,3),c(-10,10),col="lightgrey",lty=3)
lines(c(4,4),c(-10,10),col="lightgrey",lty=3)
lines(c(5,5),c(-10,10),col="lightgrey",lty=3)
lines(c(6,6),c(-10,10),col="lightgrey",lty=3)
boxplot(toto$Y~toto$X,col=c("grey","violetred1","darkolivegreen1","coral","cadetblue1","darkblue"),axes=FALSE,add=TRUE)
lines(c(-1,10),c(1,1))
box()
## Sigma
toto = data.frame(c(valid_method1$sigmaML/sigma,valid_method2$sigmaML/sigma,valid_method3$sigmaML/sigma,valid_method4$sigmaML/sigma,valid_method5$sigmaML/sigma,valid_method6$sigmaML/sigma),
                  c(seq.int(from=1, to=1, length.out=1000),seq.int(from=2, to=2, length.out=1000),
                    seq.int(from=3, to=3, length.out=1000),seq.int(from=4, to=4, length.out=1000),
                    seq.int(from=5, to=5, length.out=1000),seq.int(from=6, to=6, length.out=1000)))
names(toto) = c("Y","X")
plot(c(1,1),c(-10,10),
     type="l",col="lightgrey",lty=3, xlim=c(0.5,6.5),ylim=c(0.5,2),xlab="Method",ylab=TeX('$\\hat{\\sigma}^{ML}/\\sigma$'),axes=FALSE)
axis(2)
axis(1,at=c(1:6),label=c("1","2","3","4","5","6"))
lines(c(2,2),c(-10,10),col="lightgrey",lty=3)
lines(c(3,3),c(-10,10),col="lightgrey",lty=3)
lines(c(4,4),c(-10,10),col="lightgrey",lty=3)
lines(c(5,5),c(-10,10),col="lightgrey",lty=3)
lines(c(6,6),c(-10,10),col="lightgrey",lty=3)
boxplot(toto$Y~toto$X,col=c("grey","violetred1","darkolivegreen1","coral","cadetblue1","darkblue"),axes=FALSE,add=TRUE)
lines(c(-1,10),c(1,1))
box()
## Xi
toto = data.frame(c(valid_method1$xiML/xi,valid_method2$xiML/xi,valid_method3$xiML/xi,valid_method4$xiML/xi,valid_method5$xiML/xi,valid_method6$xiML/xi),
                  c(seq.int(from=1, to=1, length.out=1000),seq.int(from=2, to=2, length.out=1000),
                    seq.int(from=3, to=3, length.out=1000),seq.int(from=4, to=4, length.out=1000),
                    seq.int(from=5, to=5, length.out=1000),seq.int(from=6, to=6, length.out=1000)))
names(toto) = c("Y","X")
plot(c(1,1),c(-10,10),
     type="l",col="lightgrey",lty=3, xlim=c(0.5,6.5),ylim=c(0.5,2),xlab="Method",ylab=TeX('$\\hat{\\xi}^{ML}/\\xi$'),axes=FALSE)
axis(2)
axis(1,at=c(1:6),label=c("1","2","3","4","5","6"))
lines(c(2,2),c(-10,10),col="lightgrey",lty=3)
lines(c(3,3),c(-10,10),col="lightgrey",lty=3)
lines(c(4,4),c(-10,10),col="lightgrey",lty=3)
lines(c(5,5),c(-10,10),col="lightgrey",lty=3)
lines(c(6,6),c(-10,10),col="lightgrey",lty=3)
boxplot(toto$Y~toto$X,col=c("grey","violetred1","darkolivegreen1","coral","cadetblue1","darkblue"),axes=FALSE,add=TRUE)
lines(c(-1,10),c(1,1))
box()

# Table 5 : Average width of the posterior credibility interval for the 100-year quantile with the Bayesian MCMC procedure for methods 1, 2, 3 and 4.
mean(valid_method1$x95 - valid_method1$x05)
mean(valid_method2$x95 - valid_method2$x05)
mean(valid_method3$x95 - valid_method3$x05)
mean(valid_method4$x95 - valid_method4$x05)
mean(valid_method5$x95 - valid_method5$x05)

# Figure 6 : Uniformity test for the credibility intervals computed with the Bayesian MCMC procedure for methods 1, 2, 3 and 4.
hist(valid_method1$excess_proba,freq=FALSE,xlab="",ylab="",main="",ylim=c(0,3),col="lightgrey",cex.axis=1.5)
lines(c(0,1),c(1,1),lty=2,lwd=2)
lines(c(0,1),c(1.5,1.5),lty=2,col="lightgrey")
lines(c(0,1),c(2,2),lty=2,col="lightgrey")
lines(c(0,1),c(2.5,2.5),lty=2,col="lightgrey")
lines(c(0,1),c(3,3),lty=2,col="lightgrey")

# Figure 7 : 90% posterior skew surge credibility intervals based on the systematic data (grey) and on the historic data with the proposed method (red) and in the ideal case (black)
x = read.table(paste("./output/method1/",site,"_adjustments_method1.txt",sep=""))
T = x$T ; IC_low_1 = x$IC_low ; IC_up_1 = x$IC_up ; x_ML_1 = x$x_ML
x = read.table(paste("./output/method3/",site,"_adjustments_method3.txt",sep=""))
IC_low_3 = x$IC_low ; IC_up_3 = x$IC_up ; x_ML_3 = x$x_ML
## Brest 
x = read.table(paste("./output/method2/",site,"_adjustments_method2.txt",sep=""))
IC_low_2 = x$IC_low ; IC_up_2 = x$IC_up ; x_ML_2 = x$x_ML
A = length(skewSurgePOT_complete)
p_i = seq(0,0,length.out=A)
for (i in 1:A) {
  p_i[i] = i/(A + 1)*lambda
}
period = c(seq(1,1,length=length(skewSurgePOT_complete)-length(skewSurgePOT)),seq(2,2,length.out=length(skewSurgePOT)))
data_fra = data.frame(skewSurgePOT_complete,period)
data_fra = data_fra[order(data_fra$skewSurgePOT_complete,decreasing=TRUE),]
data_fra[,3] = 1/p_i
plot(T,IC_low_1,col="grey",type="l",lwd=2,
     log="x",xlim=c(1,1000),ylim=c(0.6,3.5),xlab="",ylab="",axes=F)
points(data_fra$V3[data_fra$period==1],data_fra$skewSurgePOT_complete[data_fra$period==1],col="red")
points(data_fra$V3[data_fra$period==2],data_fra$skewSurgePOT_complete[data_fra$period==2],pch=19)
lines(T,IC_low_1,lwd=2,col="grey") ; lines(T,IC_up_1,lwd=2,col="grey")
lines(T,IC_low_2,lwd=2) ; lines(T,IC_up_2,lwd=2)
lines(T,IC_low_3,lwd=2,col="red") ; lines(T,IC_up_3,lwd=2,col="red")
axis(1,at=c(1,10,100,1000),label=c("1","10","100","1000"))
axis(2,at=c(1,2,3),label=c("1","2","3"))
legend("topleft",
       legend=c("Systematic skew surges","Skew surges (1846-1952)","Method 1","Method 2","Method 3"),
       col=c("black","red","grey","black","red"),pch=c(19,1,NA,NA,NA),lty=c(NA,NA,1,1,1),lwd=c(NA,NA,2,2,2))
box()
## Dunkerque
plot(T,IC_low_1,lwd=2,col="grey",type="l",
     log="x",xlim=c(1,1000),ylim=c(0.6,3.5),xlab="",ylab="",axes=F)
# Empirical probabilities of Hirsch and Stedinger (1987)
n = length(xS) ; h_x = length(xH)
w = c(seq(0,0,length.out=n),seq(1,1,length.out=h_x)) ; xx = c(xS,xH) 
N = length(xx)
A = length(xx[xx>=u_H]) ; B = length(xx) ; C = w_H*33.33/100*lambda - h_x
p_c = A/(A + B + C) # Conditionnal probability
p_e = p_c           # Exceeding probability of the threshold
p_i = matrix(0,nrow=1,ncol=N) # Empirical probabilities
ii = 1
for (i in 1:A) {
  p_i[ii] = p_e*(i - alpha)/(A + 1 - 2*alpha)
  ii = ii + 1
}
A = length(xx[xx<u_H])
for (i in 1:A) {
  p_i[ii] = p_e + (1-p_e)*(i - alpha)/(A + 1 - 2*alpha)*lambda
  ii = ii + 1
}
TT = sort(c(1/p_i),decreasing=TRUE)
data = data.frame(xx,w) ; data = data[order(data$xx,decreasing=TRUE),]
data[,3] = TT
points(data[data$w==0,3],data[data$w==0,1],pch=19)
points(data[data$w==1,3],data[data$w==1,1],pch=1,col="red")
lines(T,IC_low_1,lwd=2,col="grey") ; lines(T,IC_up_1,lwd=2,col="grey")
lines(T,IC_low_3,lwd=2,col="red") ; lines(T,IC_up_3,lwd=2,col="red")
axis(1,at=c(1,10,100,1000),label=c("1","10","100","1000"),cex.axis=1.5)
axis(2,at=c(1,2,3),label=c("1","2","3"),cex.axis=1.5)
legend("topleft",
       legend=c("Systematic skew surges","Historical skew surges", "Method 1", "Method 3"),
       col=c("black","red","grey","red"),pch=c(19,1,NA,NA),lty=c(NA,NA,1,1),lwd=c(NA,NA,2,2),cex=1.5)
box()
## La Rochelle
plot(T,IC_low_1,lwd=2,col="grey",type="l",
     log="x",xlim=c(1,1000),ylim=c(0.6,3.5),xlab="",ylab="",axes=F)
# Empirical probabilities of Hirsch and Stedinger (1987)
n = length(skewSurgePOT) ; h_x = length(xH)
w = c(seq(0,0,length.out=n),seq(1,1,length.out=h_x)) ; xx = c(skewSurgePOT,xH) 
N = length(xx)
A = length(xx[xx>=u_H]) ; B = length(xx) ; C = w_H*8.33/100*lambda - h_x
p_c = A/(A + B + C) # Conditionnal probability
p_e = p_c           # Exceeding probability of the threshold
p_i = matrix(0,nrow=1,ncol=N) # Empirical probabilities
ii = 1
for (i in 1:A) {
  p_i[ii] = p_e*i/(A + 1)
  ii = ii + 1
}
A = length(xx[xx<u_H])
for (i in 1:A) {
  p_i[ii] = p_e + (1-p_e)*i/(A + 1)*lambda
  ii = ii + 1
}
TT = sort(c(1/p_i),decreasing=TRUE)
data = data.frame(xx,w) ; data = data[order(data$xx,decreasing=TRUE),]
data[,3] = TT
points(data[data$w==0,3],data[data$w==0,1],pch=19)
points(data[data$w==1,3],data[data$w==1,1],pch=1,col="red")
lines(T,IC_low_1,lwd=2,col="grey") ; lines(T,IC_up_1,lwd=2,col="grey")
lines(T,IC_low_3,lwd=2,col="red") ; lines(T,IC_up_3,lwd=2,col="red")
axis(1,at=c(1,10,100,1000),label=c("1","10","100","1000"),cex.main=1.5)
axis(2,at=c(1,2,3),label=c("1","2","3"),cex.main=1.5)
legend("topleft",
       legend=c("Systematic skew surges","Historical skew surges","Method 1","Method 3"),
       col=c("black","red","grey","red"),
       pch=c(19,1,NA,NA),lty=c(NA,NA,1,1),lwd=c(NA,NA,2,2),cex=1.5)
box()
## Saint Nazaire
x = read.table(paste("./output/method2/",site,"_adjustments_method2.txt",sep=""))
IC_low_2 = x$IC_low ; IC_up_2 = x$IC_up ; x_ML_2 = x$x_ML
A = length(skewSurgePOT_complete)
p_i = seq(0,0,length.out=A)
for (i in 1:A) {
  p_i[i] = i/(A + 1)*lambda
}
period = c(seq(1,1,length=length(skewSurgePOT_complete)-length(skewSurgePOT)),seq(2,2,length.out=length(skewSurgePOT)))
data_fra = data.frame(skewSurgePOT_complete,period)
data_fra = data_fra[order(data_fra$skewSurgePOT_complete,decreasing=TRUE),]
data_fra[,3] = 1/p_i
plot(T,IC_low_1,col="grey",type="l",lwd=2,
     log="x",xlim=c(1,1000),ylim=c(0.6,3.5),xlab="",ylab="",axes=F)
points(data_fra$V3[data_fra$period==1],data_fra$skewSurgePOT_complete[data_fra$period==1],col="red")
points(data_fra$V3[data_fra$period==2],data_fra$skewSurgePOT_complete[data_fra$period==2],pch=19)
lines(T,IC_low_1,lwd=2,col="grey") ; lines(T,IC_up_1,lwd=2,col="grey")
lines(T,IC_low_2,lwd=2) ; lines(T,IC_up_2,lwd=2)
lines(T,IC_low_3,lwd=2,col="red") ; lines(T,IC_up_3,lwd=2,col="red")
axis(1,at=c(1,10,100,1000),label=c("1","10","100","1000"))
axis(2,at=c(1,2,3),label=c("1","2","3"))
legend("topleft",
       legend=c("Systematic skew surges","Skew surges (1863-1956)","Method 1","Method 2","Method 3"),
       col=c("black","red","grey","black","red"),pch=c(19,1,NA,NA,NA),lty=c(NA,NA,1,1,1),lwd=c(NA,NA,2,2,2))
box()

# Table F1 : 1000-year quantile estimations obtained from the real datasets with methods 1, 2 (only for Brest and Saint Nazaire) and 3.
x1 = c(IC_low_1[1008],x_ML_1[1008],IC_up_1[1008],IC_up_1[1008] - IC_low_1[1008],(IC_up_1[1008] - IC_low_1[1008])/x_ML_1[1008]*100)
x2 = c(IC_low_2[1008],x_ML_2[1008],IC_up_2[1008],IC_up_2[1008] - IC_low_2[1008],(IC_up_2[1008] - IC_low_2[1008])/x_ML_2[1008]*100) # For Brest and Saint Nazaire
x3 = c(IC_low_3[1008],x_ML_3[1008],IC_up_3[1008],IC_up_3[1008] - IC_low_3[1008],(IC_up_3[1008] - IC_low_3[1008])/x_ML_3[1008]*100)
x = read.table(paste("./output/method3b/",site,"_adjustments_method3b.txt",sep=""))
IC_low_3b = x$IC_low ; IC_up_3b = x$IC_up ; x_ML_3b = x$x_ML
x3b = c(IC_low_3b[1008],x_ML_3b[1008],IC_up_3b[1008],IC_up_3b[1008] - IC_low_3b[1008],(IC_up_3b[1008] - IC_low_3b[1008])/x_ML_3b[1008]*100)
