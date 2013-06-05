rm(list=ls())
gene = 1
RDatafile <- paste("~/BackupProEvo/Newton/rokas7/simtree/simMax/gene",gene,".RData",sep="")
load(RDatafile)
source("~/proteinevoutk20/pkg/R/prune.R")
source("~/proteinevoutk20/pkg/R/plot_simtree.R")
ftny.vec <- apply(datanum,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
## find the starting and ending position of all edges in the tree
#plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=TRUE,plotdis=FALSE)
plot_trace(sim=sim$new[[1]],plottip=TRUE,plotftny=TRUE,plotdis=FALSE)
#################################################
#plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=FALSE,plotdis=TRUE)
plot_trace(sim=sim$new[[1]],plottip=TRUE,plotftny=FALSE,plotdis=TRUE)