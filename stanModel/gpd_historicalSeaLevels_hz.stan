functions {
  
  real gpd_lpdf (vector xS, real w_S, real u, real lambda, real sigma, real xi, real Lh) {
    // Returns the log-likelihood of the systematic skew surges and historical sea levels 

    int n;
    n = rows(xS);
    
    if (sigma<= 1e-15) {
      reject("sigma<=1e-15; found sigma =", sigma)
    }
    
    if (lambda<= 1e-15) {
      reject("lambda<=1e-15; found lambda =", lambda)
    }
    
    if (fabs(xi) > 1e-15) {
      return (n*log(lambda*w_S) - lambda*w_S - (1 + inv(xi))*sum(log(1 + (xS - u)*(xi/sigma))) - n*log(sigma) + Lh);
      
    } else {
      return (n*log(lambda*w_S) - lambda*w_S - sum((xS - u)/sigma) - n*log(sigma) + Lh);
    }
    
  }
  
  real gpd_cdf (vector x, matrix tide_d, real u, real sigma, real xi) {
    // Returns the value of the cumulative distribution function of the sea levels when their corresponding skew surges are larger than u
    
    int nn; 
    real res;
    
	  nn = rows(tide_d);
    res = 0;
    
    if (fabs(xi) > 1e-15) {

      for (j in 1:nn) {
		    if (x[j]>u) {
		      res = res + (1 - pow((1 + (x[j] - u)*xi/sigma),(- 1./xi)))*tide_d[j,2];
		    }
      }
      
    } else {
      for (j in 1:nn) {
		    if (x[j]>u) {
		      res = res + (1 - exp(- (x[j] - u)/sigma))*tide_d[j,2];
		    }
      }
    }
    
    return res;
  }
  
}

data {
  // SYSTEMATIC INFORMATION (SKEW SURGES)
  real u;                 // POT threshold of the systematic skew surges 
  int<lower=0> n;         // Length of the systematic skew surges POT sample
  vector<lower=u>[n] xS;  // Systematic skew surges POT sample
  real<lower=0> w_S;      // Systematic duration
  // EMPIRICAL DISTRIBUTION OF TIDES
  int nn;                // Number of tide intervals  
  matrix[nn,2] tide_d;   // Probability and value of tide intervals
  real p;                // Probability that a skew surge is lower than the POT threshold
  // HISTORICAL INFORMATION (SEA LEVELS)
  vector[nn] eta_H_less_tide_d;  // Historical threshold (eta_H) - tide
  real P_eta_H;                  // Probability of non exceedance of eta_H when the skew surge is lower than u
  real N;                        // Number of historical sea levels (lower and larger than eta_H) during the historical period
  int<lower=0> h_z;              // Number of historical sea levels larger than eta_H during the historical period 
  matrix[nn,h_z] xH_low;         // Historical sea levels (larger than eta_H) - tide
  matrix[nn,h_z] xH_up;          // [Historical sea levels (larger than eta_H)]*1.01 - tide
  vector[h_z] P_xH_low;          // Probability of non exceedance of xH_low
  vector[h_z] P_xH_up;           // Probability of non exceedance of xH_up
}

transformed data {
  real xmax[2];
  real MAX;
  xmax[1] = max(xH_up);
  xmax[2] = max(xS);
  MAX = max(xmax);
}

parameters {
  real<lower=0> lambda;           // Poisson parameter
  real<lower=0> sigma;            // Scale parameter of the GP distribution
  real<lower=-sigma/(MAX-u)> xi;  // Shape parameter of the GP distribution 
}

model {
  real Lh;
  real Ph;
  vector[h_z] g;
  
  // Log-probability of observing (N - h_z) sea levels not exceeding eta_H during the historical duration w_H (years)
  Ph = (N - h_z)*log(p*P_eta_H + (1 - p)*gpd_cdf(eta_H_less_tide_d, tide_d, u, sigma, xi));
  
  // Densities of the historical sea levels exceeding eta_H
  for (i in 1:h_z) {
    g[i] = ((p*P_xH_up[i] + (1 - p)*gpd_cdf(xH_up[,i], tide_d, u, sigma, xi)) - (p*P_xH_low[i] + (1 - p)*gpd_cdf(xH_low[,i], tide_d, u, sigma, xi)))/(xH_up[1,i] - xH_low[1,i]);
  }
  
  // Log-likelihood of the historical sea levels
  Lh = Ph + sum(log(g));
  
  // Log-likelihood of the systematic skew surges and historical sea levels
  xS ~ gpd(w_S, u, lambda, sigma, xi, Lh);
}
