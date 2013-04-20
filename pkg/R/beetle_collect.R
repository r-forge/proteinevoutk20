#maxdir <- "~/BackupProEvo/Newton/beetle/tree1/"
maxdir <- "~/BackupProEvo/Newton/beetle/tree2/"
l <- 8
res_max <- vector("list",length=8)

for(genect in 1:l){
  filename = paste(maxdir, "bee_",genect,".RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  else{
    load(filename)
    res_max[[genect]] <- bee.optim
  }
}

s_max <- sapply(1:8,function(x) res_max[[x]]$s)
loglik_max <- sapply(1:8,function(x) res_max[[x]]$ll$loglik)
GM_max <- sapply(1:8,function(x) res_max[[x]]$GMweights)
Q_max <- sapply(1:8,function(x) res_max[[x]]$Q)
br_max <- sapply(1:8,function(x) sum(res_max[[x]]$tree$edge.length))

#sw_max <- rbind(s_max,GM_max[2,],GM_max[3,])