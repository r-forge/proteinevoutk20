## import data on 7 tips and 8 tips
#filename: fasta file with DNA data
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
source("~/proteinevoutk20/pkg/R/getAAmodels.R")


#######################################################################
gene = 1
filename <- paste("~/proteinevoutk20/pkg/Data/Rokas/gene",gene,".fasta",sep="")
em_model <- "LG+G+F"
model <- strsplit(em_model,"+",fixed=TRUE)[[1]] #arguments in the model
p1 <- prune_new(filename,2,ROKAS_TREE)
p2 <- prune_emp(filename,2,ROKAS_TREE,em_model)
## do simulations on the regrafted branch, under both new model and the best empirical model
nsim <- 2
##simulation under both models, starting from ancestral states inferred from both models
## sim_emp_new: start from emp model result and do simulation under new model
sim_emp_emp <- vector(mode="list",length=nsim)
sim_emp_new <- vector(mode="list",length=nsim)
sim_new_emp <- vector(mode="list",length=nsim)
sim_new_new <- vector(mode="list",length=nsim)
sim_emp_emp_info <- vector(mode="list",length=nsim)
sim_emp_new_info <- vector(mode="list",length=nsim)
sim_new_emp_info <- vector(mode="list",length=nsim)
sim_new_new_info <- vector(mode="list",length=nsim)

prob_emp <- p2$prob
prob_new <- p1$prob
brlen_new <- p1$brlen
brlen_emp <- p2$brlen
s=p1$res$s
dismat = p1$res$dismat
mumat = p1$res$mumat
inv = p2$res$inv
shape=p2$res$shape
k=4
inv = 0
bf=p1$res$bfaa
#bf = NULL
opaa_p <- p1$res$ll$opaa
# Qall_p <- p1$res$Qall
index_p <- attr(p1$res$data,"index")
roots <- vector(mode="list",length=nsim)
# brlen_new = brlen_emp
# shape=1.451
for(i in 1:nsim){
  if(exists(".Random.seed"))
    rm(.Random.seed)
  root_emp <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_emp[index_p[x],]))
  root_new <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_new[index_p[x],]))
  roots[[i]]$emp <- root_emp
  roots[[i]]$new <- root_new
  sim_emp_new[[i]] <- simulation(protein=root_emp,protein_op=opaa_p[index_p],t=brlen_emp,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim_emp_emp[[i]] <- simAA(rootseq=root_emp,t=brlen_emp,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
  sim_new_new[[i]] <- simulation(protein=root_new,protein_op=opaa_p[index_p],t=brlen_new,s=s,DisMat=dismat,MuMat=mumat,bfaa=bf)
  sim_new_emp[[i]] <- simAA(rootseq=root_new,t=brlen_new,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
}

for(i in 1:nsim){
  sim_emp_new_info[[i]] <- sim.info(sim_emp_new[[i]],opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_emp_emp_info[[i]] <- sim.info(sim_emp_emp[[i]]$seq,opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_new_new_info[[i]] <- sim.info(sim_new_new[[i]],opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_new_emp_info[[i]] <- sim.info(sim_new_emp[[i]]$seq,opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
}

#sim <- simulation1(protein=root,protein_op=opaa_p[index_p],t=brlen,matall=Qall_p)
#siminfo <- sim.info(sim=sim$path,opaa=opaa_p[index_p],s=s_p,beta=beta_p,gamma=gamma_p)
brlen <- max(brlen_emp,brlen_new)
ftylim = range(c(sim_emp_emp_info[[1]]$fty,sim_emp_new_info[[1]]$fty,sim_emp_new_info[[1]]$fty,sim_new_new_info[[1]]$fty))
plot(sim_new_emp_info[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),xlim=c(0,brlen),ylim=ftylim,pch=20,xaxs="i")
for(i in 1:nsim){
  plot(sim_new_emp_info[[i]]$ftyfun,pch=20,add=TRUE)
  plot(sim_new_new_info[[i]]$ftyfun,pch=20,col="blue",add=TRUE)
  plot(sim_emp_emp_info[[i]]$ftyfun,pch=20,col="green",add=TRUE)
  plot(sim_emp_new_info[[i]]$ftyfun,pch=20,col="red",add=TRUE)
}
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
# ##Grantham distances from observed sequence
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
