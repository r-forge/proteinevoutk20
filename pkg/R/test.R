rm(list=ls())
gene = 60
RDatafile <- paste("~/BackupProEvo/Newton/rokas/simtree/gene",gene,".RData",sep="")
load(RDatafile)
source("~/proteinevoutk20/pkg/R/prune.R")
source("~/proteinevoutk20/pkg/R/plot_simtree.R")
obs.data <- conv(filename=fastafile,type="num") #observed data sequences
ftny.vec <- apply(obs.data,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
## find the starting and ending position of all edges in the tree
plot_trace(sim=sim$emp[[1]],model="emp",plottip=FALSE)
plot_trace(sim=sim$new[[1]],model="new")
