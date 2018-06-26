rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(ggplot2)

setwd("~/Git/mymsd/hapgen2_macosx_intel/")

# HAPGEN2 inputs

legs = read.table("ex.leg", header = T)
maps = read.table("ex.map", header = T)
haps = read.table("ex.haps")

plot(maps$position, maps$COMBINED_rate.cM.Mb.)
plot(maps$position, maps$Genetic_Map.cM.)

# centimorgan (c.M.) or map unit (m.u.) is a measure of recombination frequency
# it is used to imply distance along the chromesome and take into account how often recommidation occurs in a region.

# causal variants: (1) 1085679 1 1.5 2.25 (2) 2190692 0 2 4

(cv1 = which(legs$position == 1085679))
(cv1 = which(legs$position == 2190692))

table(as.matrix(haps)[cv1, ])
table(as.matrix(haps)[cv1 + 1, ])
table(as.matrix(haps)[cv1 - 1, ])

# HAPGEN2 outputs

dir("./test")

snps = read.table("test/ex.out.legend", header = T)

g.case = read.table("test/ex.out.cases.haps")
g.ctrl = read.table("test/ex.out.controls.haps")

index = 1:100 * 2
index = cbind(index - 1, index)

g.case = apply(index, 1, function(x) rowSums(g.case[, x])) %>% as.data.frame
g.ctrl = apply(index, 1, function(x) rowSums(g.ctrl[, x])) %>% as.data.frame

p.case = read.table("test/ex.out.cases.sample", header = T)
p.case = p.case[-1, ] %>% as.data.frame

p.ctrl = read.table("test/ex.out.controls.sample", header = T)
p.ctrl = p.ctrl[-1, ] %>% as.data.frame

y = c(p.ctrl$pheno, p.case$pheno) %>% as.numeric
x = cbind(g.ctrl, g.case)

# Regression

xvar = apply(x, 1, var)
x2 = x[xvar != 0, ]
snps2 = snps[xvar != 0, ]

myfit = apply(x2, 1, function(snp) {
  f = glm(y ~ snp, family = binomial)
  summary(f)$coefficients %>% as.data.frame
})

snps2$pval = sapply(myfit, function(fit) fit["snp", "Pr(>|z|)"])
snps3 = snps2[order(snps2$pval), ]
snps3$qval = snps3$pval * nrow(snps3)

plot(snps2$pos, -log10(snps2$pval))

# Hap file has 1 more row than the legend file
# HAPGEN2 phenotypes had 0 variations for the 2 leading SNPs, why?

# understand how variants LD relates to estimates

x.maf = rowSums(x) / ncol(x) / 2
x.maf[x.maf > 0.5] = 1 - x.maf[x.maf > 0.5]

geno = x[x.maf > .05, ] # use maf cutoff
snps = snps[x.maf > .05, ]

geno = apply(geno, 1, scale) %>% as.data.frame

pheno = rnorm(200, 0, 1) # random
pheno = 2 * geno[, 2] + rnorm(200, 0, 1) # causal

summary(lm(pheno ~ geno[, 2]))
summary(lm(pheno ~ geno[, 3]))
# multivariate regression?

z1 = lm(y ~ as.matrix(geno))
z1 = summary(z1)

z2 = t(geno) %*% pheno / sqrt(200)

# t(X) * Y / sqrt(N) applies to univariate analysis

z3 = apply(geno, 2, function(x) {
  f1 = summary(lm(pheno ~ x))
  f1$coefficients["x", "t value"]
})

plot(z3, c(z2))

# what about variance-covariance matrix of the effect?
# we would expect the variance-covariance matrix of the effects and the predictors are mapping into each other

library(rstan)

mystan = "

data {
  int<lower=0> N; // number of data items
  int<lower=0> K; // number of predictors
  matrix[N, K] x; // predictor matrix
  vector[N] y; // outcome vector
}
parameters {
  real alpha; // intercept
  vector[K] beta;
  real<lower=0> sigma; // coefficients for predictors
}
model {
  y ~ normal(alpha + x * beta, sigma);  // likelihood
}

"
model = stan_model(model_code = mystan)

z4 = apply(geno, 2, function(snp) {
  data = list(N = 200, K = 1, x = as.matrix(snp, nrow = 200), y = pheno)
  myfit = sampling(model, data = data, chains = 2, iter = 2e2)
  as.data.frame(myfit, pars = "beta")
})

z4 = do.call(cbind, z4)

plot(z4[, 2], z4[, 3])

# it is not obvious how to demonstrate the co-vary property of the effect estimations
# it is not obvious how to understand and use: multivariate normal null distribution of standardized effects

library(MASS)

ss = mvrnorm(n = 1e2, mu = rep(0, nrow(geno)), Sigma = var(geno))
ss = mvrnorm(n = 1e2, mu = rep(0, 2), Sigma = matrix(c(10, 3, 3, 2), nrow = 2))
