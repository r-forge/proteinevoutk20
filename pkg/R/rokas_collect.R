##lab9
opwdir <- "~/BackupProEvo/rokas_opw/"
maxdir <- "~/BackupProEvo/rokas_max/"
majdir <- "~/BackupProEvo/rokas_maj/"
##newton/ betula
# opwdir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_opw/"
# maxdir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_max/"
# majdir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_maj/"

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
opw <- sapply(1:106,function(x) res_opw[[x]]$opw)
###########################################################
opw_Nedir <- "~/BackupProEvo/rokas_opw_Ne/"
max_Nedir <- "~/BackupProEvo/rokas_max_Ne/"
maj_Nedir <- "~/BackupProEvo/rokas_maj_Ne/"

# opw_Nedir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_opw_Ne/"
# max_Nedir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_max_Ne/"
# maj_Nedir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_maj_Ne/"

res_max_Ne <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(max_Nedir,"gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_max_Ne[[genect]] <- res_op
}

s_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$s)
loglik_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$ll$loglik)
GM_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$GMweights)
Q_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$Q)
br_max_Ne <- sapply(1:106,function(x) sum(res_max_Ne[[x]]$tree$edge.length))


res_maj_Ne <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maj_Nedir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_maj_Ne[[genect]] <- res_op
}

s_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$s)
loglik_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$ll$loglik)
GM_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$GMweights)
Q_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$Q)
br_maj_Ne <- sapply(1:106,function(x) sum(res_maj_Ne[[x]]$tree$edge.length))


res_opw_Ne <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(opw_Nedir,"gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_opw_Ne[[genect]] <- res_op
}

s_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$s)
loglik_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$ll$loglik)
GM_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$GMweights)
Q_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$Q)
br_opw_Ne <- sapply(1:106,function(x) sum(res_opw_Ne[[x]]$tree$edge.length))
opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$opw)