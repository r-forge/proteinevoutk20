source("~/proteinevoutk20/pkg/R/prune.R")
source("~/proteinevoutk20/pkg/R/plot_simtree.R")
summary.sim <- vector("list",length=106)
for(gene in 1:106){
  RDatafile <- paste("~/BackupProEvo/Newton/rokas/rootMax/gene",gene,".RData",sep="")
  load(RDatafile)
  datanum <- sapply(1:length(dataphy),function(x) as.numeric(factor(tolower(data[x,]),tolower(AA))))
  datanum <- t(datanum)
  dimnames(datanum)[[1]] <- dimnames(data)[[1]]
  res_op$tree$edge.length
  index <- attr(dataphy,"index")
  opaa <- res_op$ll$opaa
  s = res_op$s
  beta = res_op$GMweights[2]
  gamma = res_op$GMweights[3]
  scale.vec <- res_op$scale.vec
  ftny.vec <- apply(datanum,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
  ftny.vec

  rootseq <- sapply(1:length(index),function(x) sample(1:20,size=1,replace=TRUE,prob=res_op$ll$ancestral[,index[x]]))
  ftny.vec <- c(ftny.vec,Ftny_protein(protein_op=opaa[index],protein=rootseq,s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma)))
  names(ftny.vec)[length(ftny.vec)] <- "root"
  sim <- simTree(tree=res_op$tree,protein_op=opaa[index],matall=res_op$Qall,rootseq=rootseq,scale.vec=scale.vec)
  trace.info <- plot_trace(sim=sim,plottip=FALSE,plotftny=FALSE,plotdis=FALSE)
  
  summary.sim[[gene]]$ftny.vec <- ftny.vec
  summary.sim[[gene]]$trace.info <- trace.info
  summary.sim[[gene]]$datanum <- datanum
  summary.sim[[gene]]$rootseq <- rootseq
  summary.sim[[gene]]$sim <- sim
}
save.image(file="summary_sim.RData",compress=TRUE)
## find the starting and ending position of all edges in the tree
#plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=TRUE,plotdis=FALSE)
# #plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=FALSE,plotdis=TRUE)
# plot_trace(sim=sim,plottip=TRUE,plotftny=FALSE,plotdis=TRUE)
#################################################
# s = 5
# fixmatall <- fixmatAll(s,DisMat=dismat)
# Qall = QAllaa1(fixmatall,mumat)
# bf <- sapply(Qall,FUN=eqmQ)
# for(i in 1:20)
#   plot(sort(bf[,i],decreasing=TRUE),type="s",ylim=c(0,1))
#################################################
