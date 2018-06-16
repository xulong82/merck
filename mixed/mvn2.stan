data {
  int<lower=1> L;
  vector[L] y;
  cov_matrix[L] Sigma;
  real prior;
}

transformed data {
  matrix[L, L] C;
  C <- cholesky_decompose(Sigma);
}

parameters {
  real mu;
  vector[L] z; 
  real<lower=machine_precision()> sigma;
  real<lower=machine_precision()> epsilon;
}

transformed parameters {
  vector[L] u;
  u <- sigma * C * z;
}

model {
  mu ~ normal(prior, 0.5);
  sigma ~ normal(prior, 0.5);
  z ~ normal(0, 1);
  y ~ normal(mu + u, epsilon);
}

generated quantities {
  vector[L] u2;
  matrix[L, L] C2;

  u2 <- sigma * C * z;
  C2 <- C;
}

