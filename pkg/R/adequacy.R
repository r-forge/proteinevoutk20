#######################################################################
# gene <- Sys.getenv("SGE_TASK_ID")
# print(gene)
# gene <- as.numeric(gene)
test = FALSE
gene = 1 #change this to run analysis on other genes
datafile <- paste("gene",gene,".RData",sep="") #RData file to load, pruning analysis
dtip = "Smik"
source("~/proteinevoutk20/pkg/R/prune.R") #resouce the files in case there are changes after *.RData file was generated
load("~/proteinevoutk20/pkg/Data/Rokas/rokas.RData")
tree <- read.nexus("~/proteinevoutk20/pkg/Data/GTR.tre")
pdffile <- paste("gene",gene,".pdf")
data = rokas[,seq(AArange[gene,1],AArange[gene,2])]
sort(unique(c(data)))
dataphy <- phyDat(data,type="AA")
dataphy
datanum <- sapply(1:length(dataphy),function(x) as.numeric(factor(tolower(data[x,]),tolower(AA))))
datanum <- t(datanum)
dimnames(datanum)[[1]] <- dimnames(data)[[1]]
obs.seq <- datanum[dtip,] #the one sequence at the deleted tip

prottest_file <- paste("~/proteinevoutk20/pkg/Result/Prottest/Rokas/rokas_",gene,"_prottest.txt",sep="")
best_emp_model <- get_best_model(filename=prottest_file)
model <- best_emp_model$model #arguments in the model
#######################################################################
## do simulations on the regrafted branch, under both new model and the best empirical model
nsim <- 10
##simulation under both models, starting from ancestral states inferred from both models
## sim_emp_new: start from emp model result and do simulation under new model
sim <- vector(mode="list")
sim$emp_emp <- vector(mode="list",length=nsim)
sim$emp_new <- vector(mode="list",length=nsim)
sim$new_emp <- vector(mode="list",length=nsim)
sim$new_new <- vector(mode="list",length=nsim)
sim_info <- sim

p1 <- prune_new(data=dataphy,dtip=dtip,tree=tree,ancestral="eqm")
p2 <- prune_emp(data=dataphy,dtip=dtip,tree=tree,model=model)
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
bf = NULL
if("F"%in% model) bf=p1$res$bfaa
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
  sim$emp_new[[i]] <- simulation(protein=root_emp,protein_op=opaa,t=brlen_emp,s=s,DisMat=dismat,MuMat=mumat)
  sim$emp_emp[[i]] <- simAA(rootseq=root_emp,t=brlen_emp,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
  sim$new_new[[i]] <- simulation(protein=root_new,protein_op=opaa,t=brlen_new,s=s,DisMat=dismat,MuMat=mumat)
  sim$new_emp[[i]] <- simAA(rootseq=root_new,t=brlen_new,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
}

## gather information from the simulated sequences
for(i in 1:nsim){
  sim_info$emp_new[[i]] <- sim.info(sim$emp_new[[i]]$path,opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$beta,gamma=p1$res$gamma)
  sim_info$emp_emp[[i]] <- sim.info(sim$emp_emp[[i]]$path,opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$beta,gamma=p1$res$gamma)
  sim_info$new_new[[i]] <- sim.info(sim$new_new[[i]]$path,opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$beta,gamma=p1$res$gamma)
  sim_info$new_emp[[i]] <- sim.info(sim$new_emp[[i]]$path,opaa=opaa,obsaa=obs.seq,
                                    s=p1$res$s,beta=p1$res$beta,gamma=p1$res$gamma)
}
########################################################################################
## plot the functionalities and the distances (from obs seq) along the simulation path
########################################################################################
pdf(pdffile)
par(mfrow=c(2,2))

##functionalities
ftyrangeall <- sapply(sim_info,ftyrange)
obs.ftny <- Ftny_protein(protein=obs.seq,protein_op=opaa,s=s,DisMat=dismat)
ftny.vec <- apply(datanum,MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)
ftny.eqm <- eqm_ftny(opaa=opaa_p,ancestral=p1$res$ll$ancestral,s=s,DisMat=dismat,index=index_p)
ftylim <- c(min(c(ftyrangeall[1,],obs.ftny,ftny.eqm)),max(c(ftyrangeall[2,],obs.ftny,ftny.eqm)))
brlen <- max(brlen_emp,brlen_new)

# sim_ftny <- vector(mode="list")
# sim_ftny$new_new <- sapply(1:nsim,function(x) tail(sim_info$new_new[[x]]$fty,1))
# sim_ftny$new_emp <- sapply(1:nsim,function(x) tail(sim_info$new_emp[[x]]$fty,1))
# sim_ftny$emp_new <- sapply(1:nsim,function(x) tail(sim_info$emp_new[[x]]$fty,1))
# sim_ftny$emp_emp <- sapply(1:nsim,function(x) tail(sim_info$emp_emp[[x]]$fty,1))
########################################################################################
plot(sim_info$new_emp[[1]]$fty~sim_info$new_emp[[1]]$t,xlab="time",ylab="functionality",
     main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_new),ylim=ftylim,
     pch=20,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  points(sim_info$new_emp[[i]]$fty~sim_info$new_emp[[i]]$t,pch=20)
  points(sim_info$new_new[[i]]$fty~sim_info$new_new[[i]]$t,pch=20,col="red")
}
abline(h=obs.ftny,col="blue") #ftny of the observed sequence at the tip
abline(h=ftny.eqm,col="green") #equilibrium ftny
## ancestral sequences estimated under empirical models
plot(sim_info$emp_emp[[1]]$fty~sim_info$emp_emp[[1]]$t,xlab="time",ylab="functionality",
     main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_emp),ylim=ftylim,
     pch=20,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  points(sim_info$emp_emp[[i]]$fty~sim_info$emp_emp[[i]]$t,pch=20)
  points(sim_info$emp_new[[i]]$fty~sim_info$emp_new[[i]]$t,pch=20,col="red")
}
abline(h=obs.ftny,col="blue")
abline(h=ftny.eqm,col="green")
## distances
disrangeall <- sapply(sim_info,disrange)
dislim <- c(min(disrangeall[1,]),max(disrangeall[2,]))
#dis.vec <- apply(datanum,MARGIN=1,FUN=pchem_d,protein2=datanum[dtip,],DisMat=dismat)
#dis.vec <- apply(dis.vec,2,mean)
plot(sim_info$new_emp[[1]]$dis~sim_info$new_emp[[1]]$t,xlab="time",ylab="match proportion",
     main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_new),ylim=dislim,
     pch=20,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  points(sim_info$new_emp[[i]]$dis~sim_info$new_emp[[i]]$t,pch=20)
  points(sim_info$new_new[[i]]$dis~sim_info$new_new[[i]]$t,pch=20,col="red")
}

plot(sim_info$emp_emp[[1]]$dis~sim_info$emp_emp[[1]]$t,xlab="time",ylab="match proportion",
     main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_emp),ylim=dislim,
     pch=20,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  points(sim_info$emp_emp[[i]]$dis~sim_info$emp_emp[[i]]$t,pch=20)
  points(sim_info$emp_new[[i]]$dis~sim_info$emp_new[[i]]$t,pch=20,col="red")
}
dev.off()
###############################################################
save.image(file=datafile,compress=TRUE)
###############################################################

