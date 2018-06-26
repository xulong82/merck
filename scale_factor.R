p = seq(.1, .9, .1)

p1 = .001

y = rbinom(1e3, 1, p1)

e = rbinom(1e3, 1, p1)
x = y
x[e == 1] = 1 - x[e == 1]

mean(y)

f1 = glm(y ~ x, family = "binomial")
f2 = glm(y ~ x, family = "gaussian")

# f1 = glm(y ~ x - 1, family = "binomial")
# f2 = glm(y ~ x - 1, family = "gaussian")

summary(f1)
summary(f2)

# mis-specification of models can cause wrong estimations on P-values

sf = mean(y) * (1 - mean(y))

ce1 = summary(f1)$coefficients
ce2 = summary(f2)$coefficients

ce2 * 1 / sf
ce1

ce = summary(f)$coefficients
ce["x", ]

f2 = glm(y ~ scale(x) - 1, family = "binomial")
summary(f2)

ce2 = summary(f2)$coefficients
ce2[1, ]
ce2[1, ] * sf

ce[1, ]

summary(x)
sd(x)
