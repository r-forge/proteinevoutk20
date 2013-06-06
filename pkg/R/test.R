rm(list=ls())
gene = 20
RDatafile <- paste("~/BackupProEvo/Newton/rokas7/simAna/Eqm/gene",gene,".RData",sep="")
load(RDatafile)
source("~/proteinevoutk20/pkg/R/prune.R")
source("~/proteinevoutk20/pkg/R/plot_simtree.R")
ftny.vec <- apply(datanum,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
## find the starting and ending position of all edges in the tree
#plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=TRUE,plotdis=FALSE)
plot_trace(sim=sim,plottip=TRUE,plotftny=TRUE,plotdis=FALSE)
#################################################
#plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=FALSE,plotdis=TRUE)
plot_trace(sim=sim,plottip=TRUE,plotftny=FALSE,plotdis=TRUE)
# s = 5
# fixmatall <- fixmatAll(s,DisMat=dismat)
# Qall = QAllaa1(fixmatall,mumat)
# bf <- sapply(Qall,FUN=eqmQ)
# for(i in 1:20)
#   plot(sort(bf[,i],decreasing=TRUE),type="s",ylim=c(0,1))
