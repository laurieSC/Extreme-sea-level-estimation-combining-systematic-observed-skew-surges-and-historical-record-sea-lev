functions {
  
  real gpd_lpdf (vector xS, real w_S, real u, real lambda, real sigma, real xi){
    // Returns the log-likelihood of the systematic skew surges POT sample
    
    int n;
    n = rows(xS);
    
    if (sigma<= 1e-15) {
      reject("sigma<=1e-15; found sigma =", sigma)
    }
    
    if (lambda<= 1e-15) {
      reject("lambda<=1e-15; found lambda =", lambda)
    }
    
    if (fabs(xi) > 1e-15) {
      return (n*log(lambda*w_S) - lambda*w_S - (1 + inv(xi))*sum(log(1 + (xS - u)*(xi/sigma))) - n*log(sigma));
      
    } else {
      return (n*log(lambda*w_S) - lambda*w_S - sum((xS - u)/sigma) - n*log(sigma));
    }
  }
  
}

data {
  real u;                 // POT threshold of the systematic skew surges 
  int<lower=0> n;         // Length of the systematic skew surges POT sample
  vector<lower=u>[n] xS;  // Systematic skew surges POT sample
  real<lower=0> w_S;      // Systematic duration
}

transformed data {
  real MAX;
  MAX = max(xS);
}

parameters {
  real<lower=0> lambda;           // Poisson parameter
  real<lower=0> sigma;            // Scale parameter of the GP distribution
  real<lower=-sigma/(MAX-u)> xi;  // Shape parameter of the GP distribution 
}

model {
  // Log-likelihood of the systematic skew surges POT sample
  xS ~ gpd(w_S, u, lambda, sigma, xi);
}
