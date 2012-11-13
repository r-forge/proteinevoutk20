maxdir <- "~/BackupProEvo/rokas_max/"
res_max <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maxdir,"gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_max[[genect]] <- res_op
}

s_max <- sapply(1:106,function(x) res_max[[x]]$s)
loglik_max <- sapply(1:106,function(x) res_max[[x]]$ll$loglik)
GM_max <- sapply(1:106,function(x) res_max[[x]]$GMweights)
Q_max <- sapply(1:106,function(x) res_max[[x]]$Q)
br_max <- sapply(1:106,function(x) sum(res_max[[x]]$tree$edge.length))

majdir <- "~/BackupProEvo/rokas_maj/"
res_maj <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(majdir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_maj[[genect]] <- res_op
}

s_maj <- sapply(1:106,function(x) res_maj[[x]]$s)
loglik_maj <- sapply(1:106,function(x) res_maj[[x]]$ll$loglik)
GM_maj <- sapply(1:106,function(x) res_maj[[x]]$GMweights)
Q_maj <- sapply(1:106,function(x) res_maj[[x]]$Q)
br_maj <- sapply(1:106,function(x) sum(res_maj[[x]]$tree$edge.length))


opwdir <- "~/BackupProEvo/rokas_opw/"
res_opw <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(opwdir,"gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_opw[[genect]] <- res_op
}

s_opw <- sapply(1:106,function(x) res_opw[[x]]$s)
loglik_opw <- sapply(1:106,function(x) res_opw[[x]]$ll$loglik)
GM_opw <- sapply(1:106,function(x) res_opw[[x]]$GMweights)
Q_opw <- sapply(1:106,function(x) res_opw[[x]]$Q)
br_opw <- sapply(1:106,function(x) sum(res_opw[[x]]$tree$edge.length))
