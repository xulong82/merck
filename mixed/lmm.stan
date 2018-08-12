data {
  int<lower=1> N;  // sample number
  vector[N] pheno; // phenotype data 
  vector[N] mymean; // mean vector of multivariate normal 
  matrix[N, N] K; // covariance of multivariate normal model 
}

parameters {
  real a; // intercept
  real<lower=1e-4> sigma_e; // variance of gaussian models
  real<lower=1e-4> sigma_g; // variance of random effect
  vector[N] u; // random effect 
}

model {
  sigma_e ~ inv_gamma(2, 1);
  sigma_g ~ inv_gamma(2, 1);
  u ~ multi_normal(mymean, sigma_g * K);

  for (n in 1:N)
    pheno[n] ~ normal(a + u[n], sigma_e);
}

