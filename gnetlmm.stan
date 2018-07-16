data {
  int<lower=1> N;  // sample number
  vector[N] x; // predictor data 
  vector[N] y; // response data 
  cholesky_factor_cov[N] L;
}

parameters {
  real a; // intercept
  real beta; // effect size of x 
  real<lower=1e-4> sigma_e; // variance of gaussian models
  real<lower=1e-4> sigma_g; // variance of random effect
  vector[N] z; // random effect primitive
}

transformed parameters {
  vector[N] u; // random effect
  u = L * z; // random effect
}


model {
  z ~ normal(0, sigma_g);
  sigma_e ~ inv_gamma(2, 1);
  sigma_g ~ inv_gamma(2, 1);

  for (n in 1:N)
    y[n] ~ normal(a + beta * x[n] + u[n], sigma_e);
}

