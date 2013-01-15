# num_sites = sapply(1:106,function(x) length(attr(ROKAS_DATA[[x]],"index")))
# ## all parameters are estimated at the same time
# opwdir_all <- "~/BackupProEvo/Newton/rokas_opall_noopw/"
# maxdir_all <- "~/BackupProEvo/Newton/rokas_opall/"
# 
# ##newton/ betula
# # opwdir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_opw/"
# # maxdir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_max/"
# # majdir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_maj/"
# 
# res_max_all <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(maxdir_all,"gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_max_all[[genect]] <- res_op
# }
# 
# s_max_all <- sapply(1:106,function(x) res_max_all[[x]]$s)
# loglik_max_all <- sapply(1:106,function(x) res_max_all[[x]]$value)
# GM_max_all <- sapply(1:106,function(x) res_max_all[[x]]$beta_gamma)
# Q_max_all <- sapply(1:106,function(x) res_max_all[[x]]$Q)
# br_max_all <- sapply(1:106,function(x) sum(res_max_all[[x]]$br))
# #################################################################################################################
# res_opw_all <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(opwdir_all,"gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_opw_all[[genect]] <- res_op
# }
# 
# s_opw_all <- sapply(1:106,function(x) res_opw_all[[x]]$s)
# loglik_opw_all <- sapply(1:106,function(x) res_opw_all[[x]]$value)
# GM_opw_all <- sapply(1:106,function(x) res_opw_all[[x]]$beta_gamma)
# Q_opw_all <- sapply(1:106,function(x) res_opw_all[[x]]$Q)
# br_opw_all <- sapply(1:106,function(x) sum(res_opw_all[[x]]$br))
# opw_all <- sapply(1:106,function(x) res_opw_all[[x]]$opw)
#################################################################################################################
opw1dir <- "~/BackupProEvo/Newton/opwFirst_short/"
opwdir <- "~/BackupProEvo/Newton/opw_short/"
maxdir <- "~/BackupProEvo/Newton/rokas_max/"
majdir <- "~/BackupProEvo/Newton/rokas_maj/"

res_max <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maxdir, "gene",genect,"_s_weight.RData",sep="")
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
sw_max <- rbind(s_max,GM_max[2,],GM_max[3,])

res_opw1 <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(opw1dir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_opw1[[genect]] <- res_op
}

s_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$s)
loglik_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$ll$loglik)
GM_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$GMweights)
Q_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$Q)
br_opw1 <- sapply(1:106,function(x) sum(res_opw1[[x]]$tree$edge.length))
opw_opw1 <- sapply(1:106, function(x) res_opw1[[x]]$opw)
sw_opw1 <- rbind(s_opw1,GM_opw1[2,],GM_opw1[3,])

res_opw <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(opwdir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename)){
    cat("load RData for gene", genect,"failed, file does not exist","\n")
    res_opw[[genect]] <- list()
  }
  else{
  load(filename)
  res_opw[[genect]] <- res_op
  }
}

s_opw <- sapply(1:106,function(x) res_opw[[x]]$s)
loglik_opw <- sapply(1:106,function(x) res_opw[[x]]$ll$loglik)
GM_opw <- sapply(1:106,function(x) res_opw[[x]]$GMweights)
Q_opw <- sapply(1:106,function(x) res_opw[[x]]$Q)
br_opw <- sapply(1:106,function(x) sum(res_opw[[x]]$tree$edge.length))
opw_opw <- sapply(1:106, function(x) res_opw[[x]]$opw)

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
sw_maj <- rbind(s_maj,GM_maj[2,],GM_maj[3,])


###########################################################
# opw_Nedir <- "~/BackupProEvo/rokas_opw_Ne/"
# max_Nedir <- "~/BackupProEvo/rokas_max_Ne/"
# maj_Nedir <- "~/BackupProEvo/rokas_maj_Ne/"
# 
# # opw_Nedir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_opw_Ne/"
# # max_Nedir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_max_Ne/"
# # maj_Nedir <- "~/proteinevoutk20/pkg/scratch/newton/rokas_maj_Ne/"
# 
# res_max_Ne <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(max_Nedir,"gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_max_Ne[[genect]] <- res_op
# }
# 
# s_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$s)
# loglik_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$ll$loglik)
# GM_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$GMweights)
# Q_max_Ne <- sapply(1:106,function(x) res_max_Ne[[x]]$Q)
# br_max_Ne <- sapply(1:106,function(x) sum(res_max_Ne[[x]]$tree$edge.length))
# 
# 
# res_maj_Ne <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(maj_Nedir, "gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_maj_Ne[[genect]] <- res_op
# }
# 
# s_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$s)
# loglik_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$ll$loglik)
# GM_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$GMweights)
# Q_maj_Ne <- sapply(1:106,function(x) res_maj_Ne[[x]]$Q)
# br_maj_Ne <- sapply(1:106,function(x) sum(res_maj_Ne[[x]]$tree$edge.length))
# 
# 
# res_opw_Ne <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(opw_Nedir,"gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_opw_Ne[[genect]] <- res_op
# }
# 
# s_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$s)
# loglik_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$ll$loglik)
# GM_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$GMweights)
# Q_opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$Q)
# br_opw_Ne <- sapply(1:106,function(x) sum(res_opw_Ne[[x]]$tree$edge.length))
# opw_Ne <- sapply(1:106,function(x) res_opw_Ne[[x]]$opw)
# 
# 
