maxbdir <- "~/BackupProEvo/Newton/beetle/tree1/"
#maxbdir <- "~/BackupProEvo/Newton/beetle/tree2/"
l <- 8
res_maxb <- vector("list",length=8)

for(genect in 1:l){
  filename = paste(maxbdir, "bee_",genect,".RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  else{
    load(filename)
    res_maxb[[genect]] <- bee.optim
  }
}

s_maxb <- sapply(1:8,function(x) res_maxb[[x]]$s)
loglik_maxb <- sapply(1:8,function(x) res_maxb[[x]]$ll$loglik)
GM_maxb <- sapply(1:8,function(x) res_maxb[[x]]$GMweights)
Q_maxb <- sapply(1:8,function(x) res_maxb[[x]]$Q)
br_maxb <- sapply(1:8,function(x) sum(res_maxb[[x]]$tree$edge.length))

#sw_maxb <- rbind(s_maxb,GM_maxb[2,],GM_maxb[3,])