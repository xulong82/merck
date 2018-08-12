# Deep-dive the colocalization analysis using coloc package

# coloc method asks whether two genetic assocation results share common genetic causal variant in a given region

# coloc method does not tell the causal relationship between the two traits, it can't separate causality from pleiotropy
# coloc value: if two associations do not share the same variant, the gene is then not the target of the trait

# https://github.com/chr1swallace/coloc/blob/master/vignettes/vignette.Rmd

rm(list = ls())

library(coloc)
library(snpStats)

setClass("simdata", representation(df1 = "data.frame", df2 = "data.frame"))
setMethod("show", signature="simdata", function(object) { cat("pair of simulated datasets, with",ncol(object@df1)-1,"SNPs and",nrow(object@df1),"samples.\n") })

nsnps = 50
nsamples = 1e3
causals = 1:2

ntotal <- nsnps * nsamples
x1 <- matrix(rbinom(ntotal,1,0.4) + rbinom(ntotal,1,0.4), ncol=nsnps)
y1 <- rnorm(nsamples, rowSums(x1[, causals]), 2)
x2 <- matrix(rbinom(ntotal,1,0.4) + rbinom(ntotal,1,0.4), ncol=nsnps)
y2 <- rnorm(nsamples, rowSums(x2[, causals]), 2)
colnames(x1) <- colnames(x2) <- paste("s", 1:nsnps, sep="")
df1 <- cbind(y = y1, x1)
df2 <- cbind(y = y2, x2)
  
data = new("simdata", df1 = as.data.frame(df1), df2 = as.data.frame(df2))

# use snpStats

Y1 <- data@df1$y
X1 <- new("SnpMatrix",as.matrix(data@df1[, -1]) + 1) # snpStat genotypes: 1/2/3
  
head(col.summary(X1))
p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1), df = 1)
  
# use lm()

Y1 <- data@df1$y
X1 <- as.matrix(data@df1[,-1])

Y1 <- data@df2$y
X1 <- as.matrix(data@df2[,-1])

test1 <- lapply(1:ncol(X1), function(i) summary(lm(Y1 ~ X1[,i]))$coefficients[2, ])

beta <- sapply(test1, function(x) x[1])
varbeta <- sapply(test1, function(x) x[2]^2)
p <- sapply(test1, function(x) x[4])

myres = data.frame(beta, varbeta, p)

W = 0.04
beta / sqrt(varbeta)
W / (varbeta + W)

coloc.abf <- finemap.abf(dataset=list(beta = beta, varbeta = varbeta, N = nrow(X1), sdY = sd(Y1), type="quant"))
head(coloc.abf)

# qts finemapping script: only needs beta and se 

myabf <- function(beta, se, W){
  V = se**2
  r = W/(V+W)
  z2 = (beta/se)**2
  abf = ((1-r)**0.5) * exp((z2/2)*r)
  return(abf)
}

myres$abf <- myabf(beta = myres$beta, se = sqrt(myres$varbeta), W = 0.04)
myres$post = myres$abf/sum(myres$abf)

# test using significant snps only

myres2 = myres[myres$p < .05, ]
myres2$abf2 <- myabf(beta = myres2$beta, se = sqrt(myres2$varbeta), W = 0.04)
myres2$post2 = myres2$abf2/sum(myres2$abf2)

# .95 credible set
myres <- myres[order(myres$post, decreasing = T), ]
myres$cumpost <- cumsum(myres$post)
credset <- myres[1:min(which(myres$cumpost > 0.95)), ]

# colocalization analysis

Y1 <- data@df1$y
Y2 <- data@df2$y
  
X1 <- as.matrix(data@df1[, -1])
X2 <- as.matrix(data@df2[, -1])
  
test1 <- lapply(1:ncol(X1), function(i) summary(lm(Y1 ~ X1[,i]))$coefficients[2, ])
test2 <- lapply(1:ncol(X2), function(i) summary(lm(Y2 ~ X2[,i]))$coefficients[2, ])
  
p1 <- sapply(test1, "[", 4)
p2 <- sapply(test2, "[", 4)
  
maf <- colMeans(X2) / 2
  
get.beta <- function(x) {
  beta <- sapply(x,"[",1)
  varbeta <- sapply(x, "[", 2)^2
  return(list(beta = beta, varbeta = varbeta))
}

b1 <- get.beta(test1)
b2 <- get.beta(test2)

dataset1 = list(pvalues = p1, N = nrow(X1), type = "quant")
dataset2 = list(pvalues = p2, N = nrow(X2), type = "quant")

dataset1$N = 1e1
dataset2$N = 1e1

# how N and MAF info were used in the coloc method?

myres <- coloc.abf(dataset1, dataset2, MAF = maf)
print(myres[[1]])
print(myres[[2]]) 

dataset1 = list(beta = b1$beta, varbeta = b1$varbeta, N = nrow(X1), sdY = sd(Y1), type = "quant") 
dataset2 = list(beta = b2$beta, varbeta = b2$varbeta, N = nrow(X2), sdY = sd(Y2), type = "quant")

myres <- coloc.abf(dataset1, dataset2, MAF = maf)
print(myres[[1]]) 
print(myres[[2]]) 
