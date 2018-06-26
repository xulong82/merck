# Christian Benner 2016 Bioinformatics
# FINEMAP: efficient variable selection using summary data from genome-wide association studies

# the likelihood function: p(y | beta, x) is proportional to N(bhat|beta, var(bhat))

x = rnorm(1e2, 0, 3)
y = 2 * x + rnorm(1e2, 0, 3)

f = lm(y ~ x - 1)
coefficients(f)

summary(f)

beta = seq(-5, 5, .1)

obj = sapply(beta, function(beta1) {
  mean1 = x * beta1
  sd1 = sd(y - mean1)
  exp(sum(dnorm(y, mean = mean1, sd = sd1, log = T)))
})

plot(beta, obj, type = "l")

bhat1 = dnorm(beta, mean = coefficients(f)["x"], sd = summary(f)$coefficients["x", "Std. Error"])
bhat2 = dnorm(coefficients(f)["x"], mean = beta, sd = summary(f)$coefficients["x", "Std. Error"])

# N(bhat|beta, var) is equivalent to N(beta|bhat, var)

all(bhat1 == bhat2)

plot(beta, bhat1)

obj.sc = obj / max(obj)
bhat.sc = bhat1 / max(bhat1)

plot(bhat.sc, obj.sc)

plot(beta, bhat.sc, type = "l", col = "blue")
lines(beta, obj.sc, col = "red")
