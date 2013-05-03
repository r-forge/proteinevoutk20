## import data on 7 tips and 8 tips
#filename: fasta file with DNA data
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
######################################################################
## analysis under new model
######################################################################
test = TRUE
prune_new <- function(filename,dtip,tree){
  data = conv(filename,type="phyDat")
  if(is.numeric(dtip)) dtip = tree$tip.label[dtip]
  tree_p = drop.tip(tree,dtip) # new tree, after trim the tip
  data_p = subset(data,subset=tree_p$tip.label) # pruned data
  data_p = phyDat(as.character(data_p),type="AA") #re-weight pruned data
  res_p <- mllm1(data_p,tree_p,s=1,beta=be,gamma=ga,Q=NU_VEC,ancestral="eqm") #ll of pruned data
  if(test){
    res_op <- res_p
    iter="10"
  }
  else{
    res_op <- optim.mllm1(res_p,optQ=T,optBranch=T,optsWeight=T,optOpw=FALSE,
                        control=list(epsilon=1e-08,hmaxit=50,htrace=1,print_level=0,maxeval="100"),ancestral="eqm")
    iter="200"
  }
  ## MLE for the parameters
  s_p = res_op$s
  beta_p=res_op$GMweights[2]
  gamma_p = res_op$GMweights[3]
  Q_p = res_op$Q
  tree_p = res_op$tree
  index_p = attr(data_p,"index")
  opaa_p = res_op$ll$opaa
  bfaa_p = res_op$bfaa
  dismat_p = res_op$dismat
  mumat_p = res_op$mumat
  fixmatall_p = res_op$fixmatall
  Qall_p <- res_op$Qall
  
  tree = ape:::reorder.phylo(tree,"pruningwise")
  ## optimize branch lengths on 8-tip tree with parameters just found, and 8-tip data
  br <- optim.br(data,tree,maxeval=iter,print_level=1,mumat=mumat_p,fixmatall=fixmatall_p,ancestral="eqm")
  br
  tree$edge.length <- br$solution #assign the tree optimized branch lengths
  br.index <- which(tree$edge[,2]==which(tree$tip.label==dtip)) #index of edge.length of the pruned branch
  brlen <- br$solution[br.index]  # length of the re-grafted branch
  tree$edge.length[br.index] <- 0 #make that branch 0
  
  nr <- attr(data_p,"nr") #number of different patterns in 7-tip data
  # for site i, find the probability (likelihood) of all 20 states at that site
  # for every site, set the observed state to be i, and then loop through all the 20 states
  state_i <- function(x){
    datai <- data_p #pruned data
    datai$add <- as.integer(rep(x,nr)) #add the pruned tip back
    names(datai)[length(datai)] <- dtip #change the name back
    ll_i <- mllm1(datai,tree,Qall=Qall_p,opaa=opaa_p,bfaa=bfaa_p,ancestral="eqm")$ll$sitelik #ll for all sites
    return(ll_i)
  }
  sitell <- sapply(1:20,state_i) #loop through all states and put results together in a matrix
  siteprob <- exp(sitell)
  return(list(brlen=brlen,tree=tree,res=res_op,prob=siteprob))
}

######################################################################
## analysis under empirical model
######################################################################
prune_emp <- function(filename,dtip,tree,model){
  data = conv(filename,type="phyDat")
  if(is.numeric(dtip)) dtip = tree$tip.label[dtip]
  tree_p = drop.tip(tree,dtip) # new tree, after trim the tip
  data_p = subset(data,subset=tree_p$tip.label) # pruned data
  data_p = phyDat(as.character(data_p),type="AA") #re-weight pruned data
  model <- strsplit(model,"+",fixed=TRUE)[[1]] #arguments in the model
  is.G <- "G" %in% model #gamma?
  is.F <- "F" %in% model #empirical frequencies?
  is.I <- "I" %in% model #invariant sites included?
  bf = NULL
  if(is.F) bf = findBf2(data_p)
  k = 1
  if(is.G) k=4
  res_p <- pml(tree=tree_p,data=data_p,bf=bf,inv=0,k=4,shape=1,model=model[1]) #ll of pruned data
  res_op <- optim.pml(res_pm,optGamma=is.G,optInv=is.I,optEdge=TRUE) #optimize parameters
  tree_p <- res_op$tree
  
  tree = ape:::reorder.phylo(ROKAS_TREE,"pruningwise")
  res <- pml(tree=tree,data=data,bf=res_op$bf,k=4,shape=res_op$shape,
               inv=res_op$inv,model=model[1]) #ll of 8-tip data
  res_op <- optim.pml(res,optEdge=TRUE) #optimize branch lengths
  ## optimize branch lengths on 8-tip tree with parameters just found, and 8-tip data
  tree <- res_op$tree
  br.index <- which(tree$edge[,2]==which(tree$tip.label==dtip)) #index of edge.length of the pruned branch
  brlen <- tree$edge.length[br.index]  # length of the re-grafted branch
  tree$edge.length[br.index] <- 0 #make that branch 0
  
  nr <- attr(data_p,"nr") #number of different patterns in 7-tip data
  # for site i, find the probability (likelihood) of all 20 states at that site
  # for every site, set the observed state to be i, and then loop through all the 20 states
  state_i <- function(x){
    datai <- data_p #pruned data
    datai$add <- as.integer(rep(x,nr)) #add the pruned tip back
    names(datai)[length(datai)] <- dtip #change the name back
    res <- pml(tree,datai,k=k,shape=res_op$shape,
               inv=res_op$inv,model=model[1]) #ll for all sites
    res <- update(res,bf=bf)
    return(res$siteLik)
  }
  sitell <- sapply(1:20,state_i) #loop through all states and put results together in a matrix
  siteprob <- exp(sitell)
  return(list(brlen=brlen,tree=tree,res=res_op,prob=siteprob))
}
# ######################################################################
gene = 1
filename <- paste("~/proteinevoutk20/pkg/Data/Rokas/gene",gene,".fasta",sep="")
p1 <- prune_new(filename,2,ROKAS_TREE)
p2 <- prune_emp(filename,2,ROKAS_TREE,"LG+G+F")
## do simulations on the regrafted branch, under both new model and the best empirical model
nsim <- 2
  ##simulation under both models, starting from ancestral states inferred from both models
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
inv = p2$res$inv
shape=p2$res$shape
k = p2$res$k
bf=p2$res$bf
opaa_p <- p1$res$ll$opaa
Qall_p <- p1$res$Qall
index_p <- attr(p1$res$data,"index")
roots <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  if(exists(".Random.seed"))
    rm(.Random.seed)
  root_emp <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_emp[index_p[x],]))
  root_new <- sapply(1:length(index_p),function(x) sample(20,1,replace=TRUE,prob=prob_new[index_p[x],]))
  roots[[i]]$emp <- root_emp
  roots[[i]]$new <- root_new
  sim_emp_new[[i]] <- simulation1(protein=root_emp,protein_op=opaa_p[index_p],t=brlen_emp,matall=Qall_p)
  sim_emp_emp[[i]] <- simAA(rootseq=root_emp,t=brlen_emp,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
  sim_new_new[[i]] <- simulation1(protein=root_new,protein_op=opaa_p[index_p],t=brlen_new,matall=Qall_p)
  sim_new_emp[[i]] <- simAA(rootseq=root_new,t=brlen_new,bf=bf,inv=inv,rate=shape,k=k,model=model[1])
}

for(i in 1:nsim){
  sim_emp_new_info[[i]] <- sim.info(sim_emp_new[[i]]$path,opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_emp_emp_info[[i]] <- sim.info(sim_emp_emp[[i]]$seq,opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_new_new_info[[i]] <- sim.info(sim_new_new[[i]]$path,opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
  sim_new_emp_info[[i]] <- sim.info(sim_new_emp[[i]]$seq,opaa=opaa_p[index_p],
                               s=p1$res$s,beta=p1$res$GMweights[2],gamma=p1$res$GMweights[3])
}
s=p1$res$s
#sim <- simulation1(protein=root,protein_op=opaa_p[index_p],t=brlen,matall=Qall_p)
#siminfo <- sim.info(sim=sim$path,opaa=opaa_p[index_p],s=s_p,beta=beta_p,gamma=gamma_p)
brlen <- max(brlen_emp,brlen_new)
plot(sim_new_emp_info[[1]]$ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),xlim=c(0,brlen),ylim=c(0,1),pch=20,xaxs="i")

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
