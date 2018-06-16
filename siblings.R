library(ggplot2)

# causal effects of genetics to trait

var = seq(1, 10, 1)
eff = seq(1, 10, 1)

out = sapply(eff, function(x) {
  
  g = rnorm(1e3, 0, 1) # genotype
  t1 = x * g + rnorm(1e3, 0, 1) # one group of the twin pairs
  t2 = x * g + rnorm(1e3, 0, 1) # the other group of the twin pairs
  
  fit1 = summary(lm(t1 ~ g))
  fit2 = summary(lm(t2 ~ g))
  
  tcor = cor(t1, t2)
  c(fit1$adj.r.squared, fit2$adj.r.squared, tcor)
})

out = data.frame(t(out))
out$var = var

plot(var, out$X1, type = "l")
lines(out$X2, col = "red")
lines(out$X3, col = "blue")

# why siblings share 50% of the genome?

# One gene: X
# Genotypes in mother: X_m1, X_m2
# Genotypes in father: X_f1, X_f2
# Genotypes in all siblings: X_m1f1, X_m2f1, X_m1f2, X_m2f2

# simulation: multiple genes of multiple siblings

ns = 1e2 # pair of siblings
ns = 1e3 # pair of siblings

yy = replicate(ns, {
  
# ng = 1 # number of genes
# ng = 1e3 # number of genes
  ng = 1e4 # number of genes
  
  y = replicate(ng, {
    m = rbinom(2, 1, 0.5) + 1
    f = rbinom(2, 1, 0.5) + 1
    
    s = paste0("m", m, "-f", f)
    s = c(sapply(s, function(x) c(gsub("-.*", "", x), gsub(".*-", "", x))))
    sum(duplicated(s)) / 2
  })
  
  sum(y) / ng
  
})

hist(yy)
