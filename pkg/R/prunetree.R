load("~/BackupProEvo/Lab9/prunetree/gene1.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
nsim <- 2
wsim <- vector(mode="list",length=nsim)
wsimWag <- vector(mode="list",length=nsim)
wroots <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  rm(.Random.seed) #make sure the random seeds reset
  ## simulation starting with WAG equilibrium frequencies
  root <- sapply(1:length(index),function(x) sample(20,1,replace=TRUE,prob=bfaa))
  wroots[[i]] <- root
  dismat = GM_cpv(GM_CPV,al,beta,gamma)
  mumat = aa_MuMat_form(res_op$Q)
  sim <- simulation(root,opaa,t=brlen,s=s,DisMat=dismat,MuMat=mumat,bfaa=bfaa) #simulation     
  wsim[[i]] <- sim.info(sim,opaa=opaa,s=s,beta=beta,gamma=gamma)
  simwag <- simAA(rootseq=root,bf=bfaa,model="LG",inv=0,rate=0.771,k=4)
  wsimWag[[i]] <- sim.info(simwag$seq,opaa=opaa,s=s,beta=beta,gamma=gamma)
  #wsimWag[[i]] <- sim.Wag(t=brlen,protein=root,bf=bfaa,opaa=opaa,s=s,dismat=dismat)
}
##################################################################################################
## analysis of the simulated data under both models
##################################################################################################
pdf("~/proteinevoutk20/pkg/Plot/prunetree/wgene1.pdf")
##functionalities based on the parameters estimated earlier
wstart.ftny.vec <-  sapply(1:nsim,function(i) head(wsim[[i]]$fty,1))
wend.ftny.vec <- sapply(1:nsim,function(i) tail(wsim[[i]]$fty,1))
plot.density(density(wend.ftny.vec),xlab="Functionality",ylab="density",main="Our model")
abline(v=Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat))
##Grantham distances from optimal aa's
wstart.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=wroots[[i]],protein2=opaa,DisMat=GM)))
wend.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(wsim[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
plot.density(density(wend.dis.vec),xlab="Avg Distance from optima",ylab="density",main="Our model")
abline(v=mean(pchem_d(datanum[6,],opaa,DisMat=dismat)))
##Grantham distances from observed sequence
wend.dis.obs.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(wsim[[x]]$sim[,1:length(index)],1)),protein2=datanum[6,],DisMat=GM)))
plot.density(density(wend.dis.obs.vec),xlab="Avg Distance from observed seq",ylab="density",main="Our model")
##################################################################################################
##WAG model
wend.ftny.vec.wag <- sapply(1:nsim,function(i) tail(wsimWag[[i]]$fty,1))
plot.density(density(wend.ftny.vec.wag),xlab="Functionality",ylab="density",main="WAG model")
abline(v=Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat))
##Grantham distances from optimal aa's
wend.dis.vec.wag <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(wsimWag[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
plot.density(density(wend.dis.vec.wag),xlab="Avg Distance from optima",ylab="density",main="WAG model")
abline(v=mean(pchem_d(datanum[6,],opaa,DisMat=dismat)))
##Grantham distances from observed sequence
wend.dis.obs.vec.wag <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(wsimWag[[x]]$sim[,1:length(index)],1)),protein2=datanum[6,],DisMat=GM)))
plot.density(density(wend.dis.obs.vec.wag),xlab="Avg Distance from observed seq",ylab="density",main="WAG model")
##################################################################################################
dev.off()
save.image(file="wgene1.RData",compress=TRUE)