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
}

transformed data {
  real xmax[2];
  real MAX;
  xmax[1] = max(eta_H_less_tide_d);
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
  
  // Log-probability of observing N sea levels not exceeding eta_H during the historical duration w_H (years)
  Ph = N*log(p*P_eta_H + (1 - p)*gpd_cdf(eta_H_less_tide_d, tide_d, u, sigma, xi));
  
  // Log-likelihood of the historical sea levels
  Lh = Ph;
  
  // Log-likelihood of the systematic skew surges and historical sea levels
  xS ~ gpd(w_S, u, lambda, sigma, xi, Lh);
}
