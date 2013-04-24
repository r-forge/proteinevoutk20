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
# opw1dir <- "~/BackupProEvo/Newton/opwFirst_short/"
# opwdir <- "~/BackupProEvo/Newton/opw_short/"

maxdir <- "~/BackupProEvo/Newton/rokas/"
#majdir <- "~/BackupProEvo/Newton/rokas_maj/"
#################################################################
res_max <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maxdir, "gene",genect,".RData",sep="")
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

maxempdir <- "~/BackupProEvo/Newton/rokas/emp_root/"
res_maxemp <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maxempdir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_maxemp[[genect]] <- res_op
}
s_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$s)
loglik_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$ll$loglik)
GM_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$GMweights)
Q_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$Q)
br_maxemp <- sapply(1:106,function(x) sum(res_maxemp[[x]]$tree$edge.length))
sw_maxemp <- rbind(s_maxemp,GM_maxemp[2,],GM_maxemp[3,])
# 
# # res_opw1 <- vector("list",length=106)
# # l <- 106
# # for(genect in 1:l){
# #   filename = paste(opw1dir, "gene",genect,"_s_weight.RData",sep="")
# #   if(!file.exists(filename))
# #     cat("load RData for gene", genect,"failed, file does not exist","\n")
# #   load(filename)
# #   res_opw1[[genect]] <- res_op
# # }
# # 
# # s_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$s)
# # loglik_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$ll$loglik)
# # GM_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$GMweights)
# # Q_opw1 <- sapply(1:106,function(x) res_opw1[[x]]$Q)
# # br_opw1 <- sapply(1:106,function(x) sum(res_opw1[[x]]$tree$edge.length))
# # opw_opw1 <- sapply(1:106, function(x) res_opw1[[x]]$opw)
# # sw_opw1 <- rbind(s_opw1,GM_opw1[2,],GM_opw1[3,])
# # 
# # res_opw <- vector("list",length=106)
# # l <- 106
# # for(genect in 1:l){
# #   filename = paste(opwdir, "gene",genect,"_s_weight.RData",sep="")
# #   if(!file.exists(filename)){
# #     cat("load RData for gene", genect,"failed, file does not exist","\n")
# #     res_opw[[genect]] <- list()
# #   }
# #   else{
# #   load(filename)
# #   res_opw[[genect]] <- res_op
# #   }
# # }
# # 
# # s_opw <- sapply(1:106,function(x) res_opw[[x]]$s)
# # loglik_opw <- sapply(1:106,function(x) res_opw[[x]]$ll$loglik)
# # GM_opw <- sapply(1:106,function(x) res_opw[[x]]$GMweights)
# # Q_opw <- sapply(1:106,function(x) res_opw[[x]]$Q)
# # br_opw <- sapply(1:106,function(x) sum(res_opw[[x]]$tree$edge.length))
# # opw_opw <- sapply(1:106, function(x) res_opw[[x]]$opw)
# 
# res_maj <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(majdir, "gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_maj[[genect]] <- res_op
# }
# 
# s_maj <- sapply(1:106,function(x) res_maj[[x]]$s)
# loglik_maj <- sapply(1:106,function(x) res_maj[[x]]$ll$loglik)
# GM_maj <- sapply(1:106,function(x) res_maj[[x]]$GMweights)
# Q_maj <- sapply(1:106,function(x) res_maj[[x]]$Q)
# br_maj <- sapply(1:106,function(x) sum(res_maj[[x]]$tree$edge.length))
# sw_maj <- rbind(s_maj,GM_maj[2,],GM_maj[3,])
# ###########################################################
# pruneDir <- "~/BackupProEvo/Lab9/prunetree/"
# diffAA <- vector(mode="numeric",length=106)
# diffAAwag <- vector(mode="numeric",length=106)
# avgDis <- vector(mode="numeric",length=106)
# avgDisWag <- vector(mode="numeric",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(pruneDir, "gene",genect,".RData",sep="")
#   if(!file.exists(filename)){
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#     pDiff <- NA
#     pDiffwag <- NA
#     end.dis.obs.vec <- NA
#     end.dis.obs.vec.wag <- NA
#   }
#   else{
#     load(filename)
#     pDiff <- sapply(1:nsim,function(x) sum(as.numeric(tail(sim[[x]]$sim[,1:length(index)],1))!=datanum[6,]))/length(index)
#     pDiffwag <- sapply(1:nsim,function(x) sum(as.numeric(tail(simWag[[x]]$sim[,1:length(index)],1))!=datanum[6,]))/length(index)
#   }
#   diffAA[genect] <- mean(pDiff)
#   diffAAwag[genect] <- mean(pDiffwag)
#   avgDis[genect] <- mean(end.dis.obs.vec)
#   avgDisWag[genect] <- mean(end.dis.obs.vec.wag)
# }
# save(s_max,diffAA,diffAAwag,avgDis,avgDisWag,file="~/Desktop/prune.RData",compress=TRUE)