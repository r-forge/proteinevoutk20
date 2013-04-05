maxdir <- "~/BackupProEvo/Newton/mammal/"
res_max <- vector("list",length=97)
l <- 97
for(genect in 1:l){
  filename = paste(maxdir, "Nemam_",genect,".RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  else{
    load(filename)
    res_max[[genect]] <- mam.optim
  }
}

s_max <- sapply(1:97,function(x) res_max[[x]]$s)
loglik_max <- sapply(1:97,function(x) res_max[[x]]$ll$loglik)
GM_max <- sapply(1:97,function(x) res_max[[x]]$GMweights)
Q_max <- sapply(1:97,function(x) res_max[[x]]$Q)
br_max <- sapply(1:97,function(x) sum(res_max[[x]]$tree$edge.length))
beta <- NULL
gamma <- NULL
for(i in 1:97){
  if(!is.null(res_max[[i]])){
    beta <- c(beta,GM_max[[i]][2])
    gamma <- c(gamma,GM_max[[i]][3])
  }
    
}
#sw_max <- rbind(s_max,GM_max[2,],GM_max[3,])