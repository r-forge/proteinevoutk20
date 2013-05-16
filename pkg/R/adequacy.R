#######################################################################
gene = 1 #change this to run analysis on other genes
datafile <- paste("~/BackupProEvo/Newton/rokas/prunetree/rootEqm/gene",gene,".RData",sep="") #RData file to load, pruning analysis
simDatafile <- paste("gene",gene,".RData",sep="")
pdffile <- paste("~/proteinevoutk20/pkg/Plot/prunetree/sim_gene",gene,".pdf",sep="") #file where the plots are saved
#pdf(pdffile)
load(datafile)
source("~/proteinevoutk20/pkg/R/main.R") #resouce the files in case there are changes after *.RData file was generated
source("~/proteinevoutk20/pkg/R/simulation.R")
obs.data <- conv(filename=fastafile,type="num") #observed data sequences
obs.seq <- obs.data["Smik",] #the one sequence at the deleted tip
model <- best_emp_model$model #arguments in the model
#######################################################################
## do simulations on the regrafted branch, under both new model and the best empirical model
nsim <- 20
##simulation under both models, starting from ancestral states inferred from both models
## sim_emp_new: start from emp model result and do simulation under new model
sim <- vector(mode="list")
sim$emp_emp <- vector(mode="list",length=nsim)
sim$emp_new <- vector(mode="list",length=nsim)
sim$new_emp <- vector(mode="list",length=nsim)
sim$new_new <- vector(mode="list",length=nsim)
sim_info <- sim

prob_emp <- p2$prob
prob_new <- p1$prob
brlen_new <- p1$brlen
brlen_emp <- p2$brlen
s=p1$res$s
dismat = p1$res$dismat
mumat = p1$res$mumat
inv = p2$res$inv
shape=p2$res$shape
if("G"%in% model) k = 4
bf=p1$res$bfaa
opaa_p <- p1$res$ll$opaa #optimal amino acids for distinct sites
index_p <- attr(p1$res$data,"index")
opaa <- opaa_p[index_p] #optimal amino acids for all sites
roots <- vector(mode="list",length=nsim)
#######################################################################
## simulations
for(i in 1:nsim){
  if(exists(".Random.seed"))
    rm(.Random.seed)
  root_emp <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_emp[x,]))
  root_new <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_new[x,]))
  roots[[i]]$emp <- root_emp
  roots[[i]]$new <- root_new
  sim$emp_new[[i]] <- simulation(protein=root_emp,protein_op=opaa,t=brlen_emp,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim$emp_emp[[i]] <- simAA(rootseq=root_emp,t=brlen_emp,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
  sim$new_new[[i]] <- simulation(protein=root_new,protein_op=opaa,t=brlen_new,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim$new_emp[[i]] <- simAA(rootseq=root_new,t=brlen_new,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
}

## gather information from the simulated sequences
for(i in 1:nsim){
  sim_info$emp_new[[i]] <- sim.info(sim$emp_new[[i]],opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_info$emp_emp[[i]] <- sim.info(sim$emp_emp[[i]]$seq,opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_info$new_new[[i]] <- sim.info(sim$new_new[[i]],opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_info$new_emp[[i]] <- sim.info(sim$new_emp[[i]]$seq,opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
}
########################################################################################
## plot the functionalities and the distances (from obs seq) along the simulation path
########################################################################################
par(mfrow=c(2,2))

##functionalities
ftyrangeall <- sapply(sim_info,ftyrange)
ftylim <- c(min(ftyrangeall[1,]),max(ftyrangeall[2,],1))
brlen <- max(brlen_emp,brlen_new)
obs.ftny <- Ftny_protein(protein=obs.data["Smik",],protein_op=opaa,s=s,DisMat=dismat)
ftny.vec <- apply(obs.data,MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)
ftny.eqm <- eqm_ftny(opaa=opaa_p,ancestral=p1$res$ll$ancestral,s=s,DisMat=dismat,index=index_p)
sim_ftny <- vector(mode="list")
sim_ftny$new_new <- sapply(1:nsim,function(x) tail(sim_info$new_new[[x]]$fty,1))
sim_ftny$new_emp <- sapply(1:nsim,function(x) tail(sim_info$new_emp[[x]]$fty,1))
sim_ftny$emp_new <- sapply(1:nsim,function(x) tail(sim_info$emp_new[[x]]$fty,1))
sim_ftny$emp_emp <- sapply(1:nsim,function(x) tail(sim_info$emp_emp[[x]]$fty,1))
########################################################################################
plot(sim_info$new_emp[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_new),ylim=ftylim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$new_emp[[i]]$ftyfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$new_new[[i]]$ftyfun,pch=20,do.points=FALSE,col="red",add=TRUE)
}
abline(h=obs.ftny,col="blue") #ftny of the observed sequence at the tip
abline(h=ftny.eqm,col="green") #equilibrium ftny
## ancestral sequences estimated under empirical models
plot(sim_info$emp_emp[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_emp),ylim=ftylim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$emp_emp[[i]]$ftyfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$emp_new[[i]]$ftyfun,pch=20,do.points=FALSE,col="red",add=TRUE)
}
abline(h=obs.ftny,col="blue")
abline(h=ftny.eqm,col="green")
## distances
disrangeall <- sapply(sim_info,disrange)
dislim <- c(min(disrangeall[1,]),max(disrangeall[2,]))

dis.vec <- apply(obs.data,MARGIN=1,FUN=pchem_d,protein2=obs.data["Smik",],DisMat=dismat)
dis.vec <- apply(dis.vec,2,mean)
plot(sim_info$new_emp[[1]]$disfun,xlab="time",ylab="distance",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_new),ylim=dislim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$new_emp[[i]]$disfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$new_new[[i]]$disfun,pch=20,do.points=FALSE,col="red",add=TRUE)
}

plot(sim_info$emp_emp[[1]]$disfun,xlab="time",ylab="distance",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_emp),ylim=dislim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$emp_emp[[i]]$disfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$emp_new[[i]]$disfun,pch=20,do.points=FALSE,col="red",add=TRUE)
}

###############################################################
## start from optimal amino acid sequence 
par(mfrow=c(1,1))
sim_t <- 1

sites <- c(0,seq(10,length(index_p),by=100))
nsites <- length(sites)
#nsites <- 2
sim_op <- vector(mode="list")
sim_op$new <- vector(mode="list",length=nsites)
sim_op$emp <- vector(mode="list",length=nsites)
sim_op_info <- sim_op

for(i in 1:nsites){
  change_ind <- sample(1:length(index_p),sites[i],replace=FALSE)
  change_aa <- sample(20,size=sites[i],replace=TRUE)
  start_aa <- opaa
  start_aa[change_ind] <- change_aa
  print(Ftny_protein(protein=start_aa,protein_op=opaa,s=s,DisMat=dismat))
  sim_op$new[[i]] <- simulation(protein=start_aa,protein_op=opaa,t=sim_t,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim_op$emp[[i]] <- simAA(rootseq=start_aa,t=sim_t,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
  sim_op_info$new[[i]] <- sim.info(sim_op$new[[i]],opaa=opaa,obsaa=obs.seq,
                                   s=s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_op_info$emp[[i]] <- sim.info(sim_op$emp[[i]]$seq,opaa=opaa,obsaa=obs.seq,
                                   s=s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
}
plot(sim_op_info$new[[1]]$ftyfun,xlab="time",ylab="functionality",
     main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,sim_t),ylim=c(0,1),
     pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE,col="red")
for(i in 1:nsites){
  plot(sim_op_info$emp[[i]]$ftyfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_op_info$new[[i]]$ftyfun,pch=20,do.points=FALSE,col="red",add=TRUE)
}
abline(h=ftny.vec,col="blue")
dev.off()
save(sim,sim_info,sim_op,sim_op_info,p1,p2,file=simDatafile)
###############################################################
## find out which one of the gene was not run, whose .RData files don't exist
# for(i in 1:106){
#   file <- paste("gene",i,".RData",sep="")
#   if(!file %in% datafiles)
#     print(i)
# }

