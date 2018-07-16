
# define phenotype skewedness

p1 = .01 
p1 = .05 
p1 = .1
p1 = .2

y = rbinom(1e3, 1, p1)
mean(y)

# define genotype skewedness (maf)

p2 = .1 
p2 = .2

e = rbinom(1e3, 1, p2)
x = y
x[e == 1] = 1 - x[e == 1]

mean(x)
mean(y)

table(x, y)

# fit w binomial / gaussian models

f1 = glm(y ~ x, family = "binomial")
f2 = glm(y ~ x, family = "gaussian")

summary(f1)
summary(f2)

# model misspecification can make wrong estimations on p-values

sf = mean(y) * (1 - mean(y))

ce1 = summary(f1)$coefficients
ce2 = summary(f2)$coefficients

ce2 * 1 / sf
ce1

# this also breaks the beta transformation formula
