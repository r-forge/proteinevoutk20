setwd("~/proteinevoutk20/pkg/scratch/newton/rokas_max")
res_max <- vector("list",length=106)
for(genect in 1:l){
  filename = paste("gene",genect,"_s_weight.RData",sep="")
  cat("load RData for gene", genect,"\n")
  load(filename)
  res_max[[genect]] <- res_op
}

s_max <- sapply(1:106,function(x) res_max[[x]]$s)
loglik_max <- sapply(1:106,function(x) res_max[[x]]$ll$loglik)
GM_max <- sapply(1:106,function(x) res_max[[x]]$GMweights)
Q_max <- sapply(1:106,function(x) res_max[[x]]$Q)
