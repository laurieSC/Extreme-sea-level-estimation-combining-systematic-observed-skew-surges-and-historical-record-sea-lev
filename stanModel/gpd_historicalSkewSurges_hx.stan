functions {
  
  real gpd_lpdf (vector xS, real w_S, real u, real lambda, real sigma, real xi, real Lh) {
    // Returns the log-likelihood of the systematic and historical skew surges
    
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
  
  real gpd_cdf (real x, real u, real sigma, real xi) {
    // Returns the value of the cumulative GP distribution function
    
    if (fabs(xi) > 1e-15) {
      return (1 - pow((1 + (x - u)*xi/sigma),(- 1./xi)));
    } else {
      return (1 - exp(-(x - u)/sigma));
    }
    
  }
  
  real gpd_df (real x, real u, real sigma, real xi) {
    // Returns the value of the density GP distribution function
    
    if (fabs(xi) > 1e-15) {
      return (pow(1 + (x - u)*xi/sigma,(- 1./xi - 1))/sigma);
    } else {
      return (exp(-(x - u)/sigma)/sigma);
    }
    
  }
  
}

data {
  // SYSTEMATIC INFORMATION (SKEW SURGES)
  real u;                  // POT threshold of the systematic skew surges 
  int<lower=0> n;          // Length of the systematic skew surges POT sample
  vector<lower=u>[n] xS;   // Systematic skew surges POT sample
  real<lower=0> w_S;       // Systematic duration
  // HISTORICAL INFORMATION (SKEW SURGES)
  real<lower=u> u_H;          // Threshold of historical skew surges
  real w_H;                   // Historical duration
  int<lower=0> h_x;           // Number of historical skew surges
  vector<lower=u_H>[h_x] xH;  // Historical skew surges
}

transformed data {
  real MAX;
  vector[2] xmax;
  xmax[1] = max(xH);
  xmax[2] = max(xS);
  MAX = max(xmax);
}

parameters {
  real<lower=0> lambda;           // Poisson parameter
  real<lower=0> sigma;            // Scale parameter of the GP distribution
  real<lower=-sigma/(MAX-u)> xi;  // Shape parameter of the GP distribution 
}

model {
  vector[h_x] f;
  real Lh;
  real Ph;
  
  // Log-probability of observing h_X skew surges exceeding u during the historical duration w_H (years)
  Ph = log(lambda*w_H) - lambda*w_H*(1 - gpd_cdf(u_H, u, sigma, xi));
  
  // Densities of the historical skew surges
  for (i in 1:h_x) {
    f[i] = gpd_df(xH[i], u, sigma, xi);
  }
  
  // Log-likelihood of the historical skew surges
  Lh = Ph + sum(log(f));

  // Log-likelihood of the systematic and historical skew surges
  xS ~ gpd(w_S, u, lambda, sigma, xi, Lh);
}
