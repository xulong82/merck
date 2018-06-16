rm(list = ls())
setwd("~/GitHub/Statistics/")

noisedeg <- read.table("noisedeg.txt")

t(means.noise<-with(noisedeg,tapply(rt,list(subj,noise),mean)))
t(means.deg<-with(noisedeg,tapply(rt,list(subj,deg),mean)))

library(lme4)
library(lattice)

# Repeated measures regression

lmList.fm1 <- lmList(rt ~ noise|subj, noisedeg)

t.test(coef(lmList.fm1)[2])

summary(m0.lmer <- lmer(rt ~ noise + (1|subj), noisedeg))

summary(lm(rt~noise,noisedeg))

print(dotplot(ranef(m0.lmer, condVar = T)))

# --- BLUP ---

library(nlme)
data(ergoStool)

(fm0 = lm(effort ~ Type - 1, ergoStool))
fm1 = lmer(effort ~ Type - 1 + (1|Subject), ergoStool)

ranef(fm1)
VarCorr(fm1)$Subject

covar.u.y = VarCorr(fm1)$Subject[1]
