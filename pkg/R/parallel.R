require(multicore)
require(foreach)
require(snow)
require(doMC)
require(doSNOW)
##number of genes in data
l <- 106
data <- list(length=l)
##data stores the data of all genes in rokas's data
system.time(
for(i in 1:l){
  data[[i]] <- conv(paste("gene",i,".fasta",sep=""))
}
)

##parallel version of the above, read in the data in 106 genes.
## doesn't save time in this case using parallelization
# print(system.time(
#   data <- mclapply(1:l,readdata <- function(x) conv(paste("gene",x,".fasta",sep="")))
#   ))

## take the number of gene as parameter, find the MLE of s for that particular gene
# ll_gene <- function(number){
#   MLE_GTR(1,0,1e5,tree,data[[number]][,1:5],m=20,al,be,ga,mumat)
# }
#system.time(mclapply(c(2:3),ll_gene, mc.preschedule=F,mc.cores=2))
##Suppose that all other parameters are the same as the default values in the MLE_GTR_GS function.
##Look for the MLE estimates of beta and gamma for all genes, and estimates for s under these values for every gene
ll_weights <- function(x){
  beta <- x[1]
  gamma <- x[2]
  #print(system.time(mle.s <- mclapply(c(1:l),ll_gene,mc.preschedule=T))) #mc.preschedule, which to use??
  ll_gene <- function(number){
    MLE_GTR(1,0,1e5,tree,data[[number]],m=20,al,beta,gamma,mumat,optim.m=2)
  }
  #browser()
  mle.s <- mclapply(c(1:l),ll_gene,mc.preschedule=T)
  
#   if(optim.m==1)
#     fvalue <- sum(unlist(lapply(1:l,function(x) as.numeric(mle.s[[x]]$fvalues)))) #sum of f values using nlminb
#   else
  fvalue <- sum(unlist(lapply(1:l,function(x) as.numeric(mle.s[[x]]$fval))))
  return(fvalue)
}

MLE_ab <- function(optim.m=1){
  if(optim.m==1)
    RES <- optimx(c(0.1,0.001),ll_weights,lower=c(0,0),upper=c(1,1),method="nlminb",hessian=T,control=list(trace=1))
  else
    RES <- bobyqa(c(0.1,0.001),ll_weights,lower=c(0,0),upper=c(1,1),control=list(iprint=3))
  return(RES)
  }
  #mclapply(c(1:2),MLE_ab,mc.preschedule=F,mc.cores=2)
