#################################################
##Look at results from optimizaiton done separately for 106 genes, where each of them has their own values for beta and gamma

gene.dir <- "~/proteinevoutk20/pkg/Result/AllGeneRes/"
##keep mle in .RData files, and discard other objects
get.mle <- function(num){
  load(paste("gene",num,".RData",sep=""))
  #keep(mle,num,sure=T)
  save(mle,file=paste("gene",num,sep=""))
}
#lapply(1:106,get.mle)

##Extract the information of MLE of parameters for each gene
##Gradient, convergence information and other info are omitted
get.par <- function(num){
  load(paste("gene",num,sep=""))
  return(mle$par$par)
}
##Do that for all genes and put them to an array
#sapply(1:106,get.par)

#################################################
load("test2")
## extract estimators for s
get.par1 <- function(num){
  res[[num]]$par
}
## extract function values using the MLE of s
get.val <- function(num){
  res[[num]]$objective
}

par <- sapply(1:106,get.par1)
fval <- sapply(1:106,get.val)
