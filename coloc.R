# Deep-dive the colocalization analysis using coloc package

# coloc method asks whether two genetic assocation results share common genetic causal variant in a given region

# coloc method does not tell the causal relationship between the two traits, it can't separate causality from pleiotropy
# positive coloc results itself would generate false positive results in identifying mediating gene of a GWAS result
# ture value of coloc is: if two associations do not share the same variant, the gene is then not the target of the trait

# https://github.com/chr1swallace/coloc/blob/master/vignettes/vignette.Rmd

library(coloc)
library(snpStats)

setClass("simdata",
         representation(df1="data.frame",df2="data.frame"))
setValidity("simdata", function(object) {
  n <- nrow(object@df1)
  if(nrow(object@df2)!=n)
    return("nrow of '@df1' should equal nrow of '@df2'")
})
setMethod("show", signature="simdata", function(object) {
  cat("pair of simulated datasets, with",ncol(object@df1)-1,"SNPs and",nrow(object@df1),"samples.\n")
})

sim.data <- function(nsnps=50,nsamples=200,causals=1:2,nsim=1) {
  cat("Generate",nsim,"small sets of data\n")
  ntotal <- nsnps * nsamples * nsim
  X1 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y1 <- rnorm(nsamples,rowSums(X1[,causals]),2)
  X2 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y2 <- rnorm(nsamples,rowSums(X2[,causals]),2)
  colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
  df1 <- cbind(Y=Y1,X1)
  df2 <- cbind(Y=Y2,X2)
  if(nsim==1) {
    return(new("simdata",
               df1=as.data.frame(df1),
               df2=as.data.frame(df2)))
  } else {
    index <- split(1:(nsamples * nsim), rep(1:nsim, nsamples))
    objects <- lapply(index, function(i) new("simdata", df1=as.data.frame(df1[i,]),
                                             df2=as.data.frame(df2[i,])))
    return(objects)
  }
}

data <- sim.data(nsamples=1000,nsim=1)

Y1 <- data@df1$Y
X1 <- new("SnpMatrix",as.matrix(data@df1[,-1]))
p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
maf <- col.summary(X1)[,"MAF"]
