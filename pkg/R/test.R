## import data on 7 tips and 8 tips
#filename: fasta file with DNA data
# source("~/proteinevoutk20/pkg/R/main.R")
# source("~/proteinevoutk20/pkg/R/simulation.R")
# source("~/proteinevoutk20/pkg/R/getAAmodels.R")
#######################################################################
gene = 1
datafile <- paste("~/BackupProEvo/Newton/rokas/prunetree/rootEmp/gene",gene,".RData",sep="")
load(datafile)
# source("~/proteinevoutk20/pkg/R/prune.R")
# p2 <- prune_emp(fastafile,"Sklu",ROKAS_TREE,best_emp_model$model)
obs.data <- conv(filename=fastafile,type="num")
obs.seq <- obs.data["Scas",]
model <- best_emp_model$model #arguments in the model

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
opaa_p <- p1$res$ll$opaa
index_p <- attr(p1$res$data,"index")
roots <- vector(mode="list",length=nsim)

for(i in 1:nsim){
  if(exists(".Random.seed"))
    rm(.Random.seed)
  root_emp <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_emp[index_p[x],]))
  root_new <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_new[index_p[x],]))
  roots[[i]]$emp <- root_emp
  roots[[i]]$new <- root_new
  sim$emp_new[[i]] <- simulation(protein=root_emp,protein_op=opaa_p[index_p],t=brlen_emp,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim$emp_emp[[i]] <- simAA(rootseq=root_emp,t=brlen_emp,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
  sim$new_new[[i]] <- simulation(protein=root_new,protein_op=opaa_p[index_p],t=brlen_new,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim$new_emp[[i]] <- simAA(rootseq=root_new,t=brlen_new,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
}
source("~/proteinevoutk20/pkg/R/simulation.R")
for(i in 1:nsim){
  sim_info$emp_new[[i]] <- sim.info(sim$emp_new[[i]],opaa=opaa_p[index_p],obsaa=obs.seq,
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_info$emp_emp[[i]] <- sim.info(sim$emp_emp[[i]]$seq,opaa=opaa_p[index_p],obsaa=obs.seq,
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_info$new_new[[i]] <- sim.info(sim$new_new[[i]],opaa=opaa_p[index_p],obsaa=obs.seq,
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_info$new_emp[[i]] <- sim.info(sim$new_emp[[i]]$seq,opaa=opaa_p[index_p],obsaa=obs.seq,
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
}



#sim <- simulation1(protein=root,protein_op=opaa_p[index_p],t=brlen,matall=Qall_p)
#siminfo <- sim.info(sim=sim$path,opaa=opaa_p[index_p],s=s_p,beta=beta_p,gamma=gamma_p)

ftyrangeall <- sapply(sim_info,ftyrange)
ftylim <- c(min(ftyrangeall[1,]),max(ftyrangeall[2,],1))
brlen <- max(brlen_emp,brlen_new)
obs.ftny <- Ftny_protein(protein=obs.data["Scas",],protein_op=opaa_p[index_p],s=s,DisMat=dismat)
ftny.vec <- apply(obs.data,MARGIN=1,FUN=Ftny_protein,protein_op=opaa_p[index_p],s=s,DisMat=dismat)
plot(sim_info$new_emp[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_new),ylim=ftylim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$new_emp[[i]]$ftyfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$new_new[[i]]$ftyfun,pch=20,do.points=FALSE,col="red",add=TRUE)
  #plot(sim_info$op_new[[i]]$ftyfun,pch=20,do.points=FALSE,col="green",add=TRUE)
}
abline(h=obs.ftny,col="blue")
plot(sim_info$emp_emp[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_emp),ylim=ftylim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$emp_emp[[i]]$ftyfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$emp_new[[i]]$ftyfun,pch=20,do.points=FALSE,col="red",add=TRUE)
  #plot(sim_info$op_emp[[i]]$ftyfun,pch=20,do.points=FALSE,col="green",add=TRUE)
}
abline(h=obs.ftny,col="blue")

disrangeall <- sapply(sim_info,disrange)
dislim <- c(min(disrangeall[1,]),max(disrangeall[2,]))

dis.vec <- apply(obs.data,MARGIN=1,FUN=pchem_d,protein2=obs.data["Scas",],DisMat=dismat)
dis.vec <- apply(dis.vec,2,mean)
plot(sim_info$new_emp[[1]]$disfun,xlab="time",ylab="distance",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_new),ylim=dislim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$new_emp[[i]]$disfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$new_new[[i]]$disfun,pch=20,do.points=FALSE,col="red",add=TRUE)
  #plot(sim_info$op_new[[i]]$ftyfun,pch=20,do.points=FALSE,col="green",add=TRUE)
}

plot(sim_info$emp_emp[[1]]$disfun,xlab="time",ylab="distance",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,brlen_emp),ylim=dislim,pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
for(i in 1:nsim){
  plot(sim_info$emp_emp[[i]]$disfun,pch=20,do.points=FALSE,add=TRUE)
  plot(sim_info$emp_new[[i]]$disfun,pch=20,do.points=FALSE,col="red",add=TRUE)
  #plot(sim_info$op_emp[[i]]$ftyfun,pch=20,do.points=FALSE,col="green",add=TRUE)
}


###############################################################
## start from optimal amino acid sequence 
# nsim <- 1
# sim_op_new <- vector(mode="list",length=nsim)
# sim_op_emp <- sim_op_new
# sim_op_new_info <- sim_op_new
# sim_op_emp_info <- sim_op_emp
# for(i in 1:nsim){
#   sim_op_new[[i]] <- simulation(protein=opaa_p[index_p],protein_op=opaa_p[index_p],t=10,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
#   sim_op_emp[[i]] <-simAA(rootseq=opaa_p[index_p],t=10,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
#   sim_op_new_info[[i]] <- sim.info(sim_op_new[[i]],opaa=opaa_p[index_p],
#                                    s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
#   sim_op_emp_info[[i]] <- sim.info(sim_op_emp[[i]]$seq,opaa=opaa_p[index_p],
#                                    s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
# }
# plot(sim_op_emp_info[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("gene ",gene,", s=",round(s,2),sep=""),xlim=c(0,10),ylim=c(0,1),pch=20,do.points=FALSE,xaxs="i",frame.plot=FALSE)
# for(i in 1:nsim){
#   plot(sim_op_emp_info[[i]]$ftyfun,pch=20,do.points=FALSE,add=TRUE)
#   plot(sim_op_new_info[[i]]$ftyfun,pch=20,do.points=FALSE,col="red",add=TRUE)
# }
# abline(h=obs.ftny,col="blue")
###############################################################
## find out which one of the gene was not run, whose .RData files don't exist
# for(i in 1:106){
#   file <- paste("gene",i,".RData",sep="")
#   if(!file %in% datafiles)
#     print(i)
# }
###############################################################  

# plot(sim_info$new_emp[[1]]$disfun,xlab="time",ylab="distance",main=paste("distance, s=",round(s,3),sep=""),xlim=c(0,brlen),pch=20,do.points=FALSE,xaxs="i",yaxs="i",frame.plot=FALSE)
# for(i in 1:nsim){
#   plot(sim_info$new_emp[[i]]$disfun,pch=20,do.points=FALSE,add=TRUE)
#   plot(sim_info$new_new[[i]]$disfun,pch=20,do.points=FALSE,col="blue",add=TRUE)
#   plot(sim_info$emp_emp[[i]]$disfun,pch=20,do.points=FALSE,col="green",add=TRUE)
#   plot(sim_info$emp_new[[i]]$disfun,pch=20,do.points=FALSE,col="red",add=TRUE)
# }
# obs.dis <- mean(pchem_d(obs.data["Sklu",],opaa_p[index_p],DisMat=dismat))
# abline(h=obs.dis,col="purple")
# ##################################################################################################
# ## analysis of the simulated data under both models
# ##################################################################################################
# pdf("~/proteinevoutk20/pkg/Plot/prunetree/gene88.pdf")
# ##functionalities based on the parameters estimated earlier
# start.ftny.vec <-  sapply(1:nsim,function(i) head(sim[[i]]$fty,1))
# end.ftny.vec <- sapply(1:nsim,function(i) tail(sim[[i]]$fty,1))
# plot.density(density(end.ftny.vec),xlab="Functionality",ylab="density",main="Our model")
# abline(v=Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat))
# ##Grantham distances from optimal aa's
# start.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=roots[[i]],protein2=opaa,DisMat=GM)))
# end.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(sim[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
# plot.density(density(end.dis.vec),xlab="Avg Distance from optima",ylab="density",main="Our model")
# abline(v=mean(pchem_d(datanum[6,],opaa,DisMat=dismat)))
##Grantham distances from observed sequence
# end.dis.obs.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(sim[[x]]$sim[,1:length(index)],1)),protein2=datanum[6,],DisMat=GM)))
# plot.density(density(end.dis.obs.vec),xlab="Avg Distance from observed seq",ylab="density",main="Our model")
# ## number of aa differences between the predictions and the observed sequence:
# diffAA <- sum(as.numeric(tail(sim[[1]]$sim[,1:length(index)],1))!=datanum[6,])
# diffAAp <- diffAA/length(index)
# ##################################################################################################
# ##WAG model
# end.ftny.vec.emp <- sapply(1:nsim,function(i) tail(simEmp[[i]]$fty,1))
# plot.density(density(end.ftny.vec.emp),xlab="Functionality",ylab="density",main="WAG model")
# abline(v=Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat))
# ##Grantham distances from optimal aa's
# end.dis.vec.emp <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(simEmp[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
# plot.density(density(end.dis.vec.emp),xlab="Avg Distance from optima",ylab="density",main="WAG model")
# abline(v=mean(pchem_d(datanum[6,],opaa,DisMat=dismat)))
# ##Grantham distances from observed sequence
# end.dis.obs.vec.emp <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(simEmp[[x]]$sim[,1:length(index)],1)),protein2=datanum[6,],DisMat=GM)))
# plot.density(density(end.dis.obs.vec.emp),xlab="Avg Distance from observed seq",ylab="density",main="WAG model")
# ##################################################################################################
# dev.off()
# 
# # plot.grid <- function(index){
# #   xyz = grid.gen(opw_mean,index,res_op=res_op,gridnum=20)
# #   akima.xyz = interp(xyz$x,xyz$y,xyz$z)
# #   image(akima.xyz)
# #   contour(akima.xyz,add=TRUE)
# #   points(xyz,pch=3,col="blue")
# # }
# 
# # beta=be
# # gamma=ga
# # s = 2
# # Q = NU_VEC
# # dismat = GM_cpv(GM_CPV,al,beta,gamma)
# # fixmatall <- fixmatAll(s,DisMat=dismat)
# # mumat = aa_MuMat_form(Q)
# # Qall = QAllaa1(fixmatall,mumat)
# 
