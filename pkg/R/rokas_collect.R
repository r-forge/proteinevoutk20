## insert newrow to existing data frame at rth place
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
# output clean model comparison for rokas data, gene by gene
model_comp <- function(i){
  protTest_file <- paste("~/proteinevoutk20/pkg/Result/Prottest/Rokas/rokas",i,".txt",sep="")
  maxroot_file <- paste("~/BackupProEvo/Newton/rokas/rootMax/gene",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/rokas/rootEmp/gene",i,".RData",sep="")
  oproot_file <- paste("~/BackupProEvo/Newton/rokas/rootOp/gene",i,".RData",sep="")
  eqmroot_file <- paste("~/BackupProEvo/Newton/rokas/rootEqm/gene",i,".RData",sep="")
  #### protTest result on empirical models
  r1 <- read.table(protTest_file,header=TRUE,as.is=TRUE)
  names(r1)[5] <- "neg.LnL"
  r1$para <- round((r1$AIC - 2*r1$neg.LnL)/2)
  r1$AICw <- NULL
  #### root states maximized
  load(maxroot_file)
  df_bs <- length(res_op$tree$edge.length)+3+5+19 #number of free parameters, without counting op sites
  nr <- attr(res_op$data,"nr")
  r1 <- insertRow(r1,list("New+max+maxroot",0,0,round(-res_op$ll$loglik,2),df_bs+2*nr),1)
  #### root states with empirical frequencies
  load(emproot_file)
  r1 <- insertRow(r1,list("New+max+emproot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  #### root states with equilibrium frequencies
  load(eqmroot_file)
  r1 <- insertRow(r1,list("New+max+eqmroot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  #### root states with empirical frequencies
  load(oproot_file)
  r1 <- insertRow(r1,list("New+max+oproot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  r1$AIC <- 2*(r1$neg.LnL+r1$para)
  r1$deltaAIC <- r1$AIC - min(r1$AIC)
  return(r1)
}
getres <- function(i){
  maxroot_file <- paste("~/BackupProEvo/Newton/rokas/rootMax/gene",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/rokas/rootEmp/gene",i,".RData",sep="")
  oproot_file <- paste("~/BackupProEvo/Newton/rokas/rootOp/gene",i,".RData",sep="")
  eqmroot_file <- paste("~/BackupProEvo/Newton/rokas/rootEqm/gene",i,".RData",sep="")
  res <- vector(mode="list",length=4)
  load(maxroot_file)
  res[[1]] <- res_op
  load(emproot_file)
  res[[2]] <- res_op
  load(oproot_file)
  res[[3]] <- res_op
  load(eqmroot_file)
  res[[4]] <- res_op
  return(res)
}
op.comp <- function(res){
  ind <- combinations(n=4,r=2)
  opaa <- sapply(1:4,function(x) res[[x]]$ll$opaa)
  comp <- NULL
  for(i in 1:6){
    comp <- c(comp,sum(opaa[,ind[i,1]]==opaa[,ind[i,2]]))
  }
  return(list(opaa=opaa,comp=comp))
}

#################################################################
maxdir <- "~/BackupProEvo/Newton/rokas/rootEqm/"
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

# maxempdir <- "~/BackupProEvo/Newton/rokas/emp_root/"
# res_maxemp <- vector("list",length=106)
# l <- 106
# for(genect in 1:l){
#   filename = paste(maxempdir, "gene",genect,"_s_weight.RData",sep="")
#   if(!file.exists(filename))
#     cat("load RData for gene", genect,"failed, file does not exist","\n")
#   load(filename)
#   res_maxemp[[genect]] <- res_op
# }
# s_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$s)
# loglik_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$ll$loglik)
# GM_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$GMweights)
# Q_maxemp <- sapply(1:106,function(x) res_maxemp[[x]]$Q)
# br_maxemp <- sapply(1:106,function(x) sum(res_maxemp[[x]]$tree$edge.length))
# sw_maxemp <- rbind(s_maxemp,GM_maxemp[2,],GM_maxemp[3,])
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