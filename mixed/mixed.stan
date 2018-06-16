data {
  int<lower=1> C; // number of groups
  int<lower=0> N; // number of data points
  int<lower=1> D; // input dimension
  int<lower=1,upper=C> l[N]; // group labels
  matrix[N,D] x;
  vector[N] y;
}
parameters {
  vector[D] w[C]; // regression weights
  real<lower=0> sigma; // noise stdev
  vector[D] mu_w; // hyperparameters
  vector<lower=0>[D] sigma_w; // hyperparameters
}
model {
  vector[N] y_mu;
  
  mu_w ~ normal(0, 10);
  sigma_w ~ cauchy(0, 5);

  for (c in 1:C) {
    w[c] ~ normal(mu_w, sigma_w);
  }
 
  sigma ~ cauchy(0, 5);
  
  for (n in 1:N) {
    y_mu[n] = x[n, ] * w[l[n]];
  }
  y ~ normal(y_mu, sigma);
}

