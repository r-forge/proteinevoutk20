load("~/BackupProEvo/Lab9/prunetree/gene83.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
nsim <- 20
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
  simwag <- simAA(rootseq=root,bf=bfaa,model="mtArt",inv=0,rate=0.621,k=4)
  wsimWag[[i]] <- sim.info(simwag$seq,opaa=opaa,s=s,beta=beta,gamma=gamma)
  #wsimWag[[i]] <- sim.Wag(t=brlen,protein=root,bf=bfaa,opaa=opaa,s=s,dismat=dismat)
}
##################################################################################################
## analysis of the simulated data under both models
##################################################################################################
#pdf("~/proteinevoutk20/pkg/Plot/prunetree/wgene1.pdf")
##functionalities based on the parameters estimated earlier
par(mfrow=c(2,3))
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
#dev.off()
#save.image(file="wgene1.RData",compress=TRUE)
par(mfrow=c(1,2))
## plot the functionalities of protein sequence on the simulation path
wsim.fty.range <- range(unlist(sapply(1:nsim,function(x) wsim[[x]]$fty)))
wsimWag.fty.range <- range(unlist(sapply(1:nsim,function(x) wsimWag[[x]]$fty)))
fty.range <- range(c(unlist(sapply(1:nsim,function(x) wsim[[x]]$fty)),unlist(sapply(1:nsim,function(x) wsimWag[[x]]$fty))))
plot(wsim[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,brlen),xaxs="i",ylim=fty.range)
for(i in 2:nsim){
  plot(wsim[[i]]$ftyfun,pch=20,xlim=c(0,brlen),xaxs="i",ylim=fty.range,add=TRUE)
}
for(i in 1:nsim){
  plot(wsimWag[[i]]$ftyfun,pch=20,xlim=c(0,brlen),xaxs="i",ylim=fty.range,col="red",add=TRUE)
}
## plot the distances of protein sequence from optimal protein on the simulation path
wsim.dis.range <- range(unlist(sapply(1:nsim,function(x) wsim[[x]]$dis)))
wsimWag.dis.range <- range(unlist(sapply(1:nsim,function(x) wsimWag[[x]]$dis)))
dis.range <- range(c(unlist(sapply(1:nsim,function(x) wsim[[x]]$dis)),unlist(sapply(1:nsim,function(x) wsimWag[[x]]$dis))))
plot(wsim[[1]]$disfun,xlab="time",ylab="distance",main=paste("distance, s=",round(s,3),sep=""),pch=20,xlim=c(0,brlen),xaxs="i",ylim=dis.range)
for(i in 2:nsim){
  plot(wsim[[i]]$disfun,pch=20,xlim=c(0,brlen),xaxs="i",ylim=dis.range,add=TRUE)
}
for(i in 1:nsim){
  plot(wsimWag[[i]]$disfun,pch=20,xlim=c(0,brlen),xaxs="i",ylim=dis.range,col="red",add=TRUE)
}

