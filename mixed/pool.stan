data {
  int<lower=0> N; // number of data points
  int<lower=1> D; // number of regressors
  matrix[N,D] x;  // matrix of regressors 
  vector[N] y;    // target outputs
}
parameters {
  vector[D] w; // regression weights
  real<lower=0> sigma; // noise stdev
}
model {
  sigma ~ cauchy(0, 5); // half-Cauchy prior
  w ~ normal(0, 10);

  y ~ normal(x * w, sigma);
}

