##Look at results from optimizaiton done separately for 106 genes, where each of them has their own values for beta and gamma
gene.dir <- "~/proteinevoutk20/pkg/Result/AllGeneRes/"


get.mle <- function(num){
  load(paste("gene",num,".RData",sep=""))
  #keep(mle,num,sure=T)
  save(mle,file=paste("gene",num,sep=""))
}
lapply(1:106,get.mle)

##Extract the information of MLE of parameters for each gene
get.par <- function(num){
  load(paste("gene",num,sep=""))
  return(mle$par$par)
}
##Do that for all genes and put them to an array
sapply(1:106,get.par)

get.val <- function(x,y){
  file <- paste("bg.",x,".",y,sep="")
  load(file)
##  par <- sapply(1:106, function(x) res[[x]]$par)
  val <- sapply(1:106, function(x) res[[x]]$objective)
  sum(val)
}
