load("~/proteinevoutk20/pkg/scratch/lab9/gene7p/gene7p.RData")
source("~/proteinevoutk20/pkg/R/main.R")
## the objects for the pruned data/tree has suffix "p" indicating that it's for the pruned
## use the paraemters found with 7 tips, find the optimal branch lengths with 8 tips
data = res_op$data #data with 7 tips
s = res_op$s
beta= res_op$GMweights[2]
gamma = res_op$GMweights[3]
Q = res_op$Q
treep = res_op$tree
index = attr(res_op$data,"index")
opaa = res_op$ll$opaa[index] 
bfaa= res_op$bfaa
res_op <- mllm1(data=data,tree=treep,s=s,beta=beta,gamma=gamma,Q=Q,bfaa=bfaa)
gene7 <- conv("~/proteinevoutk20/pkg/Data/gene7.fasta",type="phyDat")
gene7p <- conv("~/proteinevoutk20/pkg/scratch/lab9/gene7p/gene7p.fasta",type="phyDat")
dataPnum <- matrix(unlist(gene7p),nrow=7,byrow=TRUE)

tree <- ROKAS_TREE
tree <- ape:::reorder.phylo(tree,"p")
load("~/proteinevoutk20/pkg/scratch/lab9/gene7p/gene7_8tips_br.RData")
tree$edge.length <- br$solution
tree$edge.length[10] <- 0 # set the edge leading to Scas to 0
brlen <- br$solution[10] 

nr <- attr(gene7p,"nr") #number of different patterns in the data with 7 tips
opaap <- res_op$ll$opaa # optimal amino acids at all different site patterns, instead of all sites
Qall <- res_op$Qall
#for each site pattern:

sitei <- function(i){
  datai <- dataPnum[,i] #ith distinct pattern in the truncated(pruned) data
  addj <- function(j){
    datai_j <- append(datai,values=j,after=5) #Add Scas as a tip but with external branch length equal to 0
    return(ll_site(tree,datai_j,optimal=opaap[i],s=s,Q=Q,alpha=al,beta=beta,gamma=gamma,bf=bfaa)) #loglikelihood at this site with Scas has state j
  }
  sitelike <- sapply(1:20,addj) # all 20 possible states at site i
  return(sitelike)
}

siteprob <- mclapply(1:nr,sitei)
siteprob <- simplify2array(siteprob,higher=FALSE)
siteprop <- apply(siteprob,2,function(x) x/sum(x)) #normalize each site prob so that the sum is 1, proportion to be used in simulation
save(siteprop,siteprob,file="proteinevoutk20/pkg/scratch/lab9/gene7p/siteprob.RData",compress=TRUE)

siteprop <- siteprop[,index]
root<- sapply(1:length(index),function(x) sample(20,1,replace=TRUE,prob=siteprop[,x]))
# in Amino Acid: 
# sample <- AA[sample]

nsim <- 100
sim <- vector(mode="list",length=nsim)
sim_start <- plot.sim(s=s,t=brlen,root=root,opaa=opaa,beta=beta,gamma=gamma,bfaa=bfaa,dist=FALSE,add=F)
roots <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  root<- sapply(1:length(index),function(x) sample(20,1,replace=TRUE,prob=siteprop[,x]))
  roots[[i]] <- root
  sim[[i]] <- plot.sim(s=s,t=brlen,root=root,opaa=opaa,beta=beta,gamma=gamma,bfaa=bfaa,dist=FALSE,add=T)
}
##functionalities based on the parameters estimated earlier
start.ftny.vec <-  sapply(1:nsim,function(i) head(sim[[i]]$fty,1))
end.ftny.vec <- sapply(1:nsim,function(i) tail(sim[[i]]$fty,1))
plot.density(density(end.ftny.vec),xlab="Functionality",ylab="density",main="Our model")
abline(v=Ftny_protein(protein=gene7num[6,],protein_op=opaa,s=s,DisMat=dismat))
##Grantham distances from optimal aa's
end.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(sim[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
plot.density(density(end.dis.vec),xlab="Avg Distance from optima",ylab="density",main="Our model")
##Grantham distances from observed sequence
end.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(sim[[x]]$sim[,1:length(index)],1)),protein2=gene7num[6,],DisMat=GM)))
plot.density(density(end.dis.vec),xlab="Avg Distance from observed seq",ylab="density",main="Our model")

simWag <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  simWag[[i]] <- plotWag(protein=roots[[i]],t=brlen,bf=bfaa,add=TRUE)
}
end.ftny.vec.wag <- sapply(1:nsim,function(i) tail(simWag[[i]]$fty,1))
plot.density(density(end.ftny.vec.wag),xlab="Functionality",ylab="density",main="WAG model")
abline(v=Ftny_protein(protein=gene7num[6,],protein_op=opaa,s=s,DisMat=dismat))
end.dis.vec.wag <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(simWag[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
plot.density(density(end.dis.vec.wag),xlab="Avg Distance from optima",ylab="density",main="WAG model")
##Grantham distances from observed sequence
end.dis.vec.wag <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(simWag[[x]]$sim[,1:length(index)],1)),protein2=gene7num[6,],DisMat=GM)))
plot.density(density(end.dis.vec.wag),xlab="Avg Distance from observed seq",ylab="density",main="WAG model")


for(i in 1:9){
  plot(simWag[[i]]$ftyfun,xlab="time",ylab="functionality",main="functionality",pch=20,xlim=c(0,brlen),xaxs="i")
}

for(i in 1:10){
  plot(sim[[i]]$ftyfun,xlab="time",ylab="functionality",main="functionality",pch=20,xlim=c(0,brlen),xaxs="i")
}
## simulation based on WAG model, and the plots of functionality along the path using the simulations
plotWag <- function(t=brlen,protein,bf=bfaa,add=FALSE){
  sim <- simulationQ(protein=protein,t=t,Q=WagMat,bf=bf)
  l <- dim(sim)[2] - 2
  fty <- apply(sim[,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
  dis <- apply(sim[,1:l],MARGIN=1,FUN=pchem_d,protein2=opaa,DisMat=dismat)#distance from optimal amino acids
  dis <- apply(dis,2,sum)#sum of distances for all sites
  ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
  disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
  #plot(ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
  return(list(sim=sim,fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
}

plot.sim <- function(s=1,t=10,root=root,opaa=opaa,beta,gamma,bfaa=bfaa,func=TRUE,dist=TRUE, add=FALSE){
  dismat = GM_cpv(GM_CPV,al,beta,gamma)
  mumat = aa_MuMat_form(res_op$Q)
  sim <- simulation(root,opaa,t=t,s=s,DisMat=dismat,MuMat=mumat,bfaa=bfaa) #simulation
  l <- dim(sim)[2] - 2
  fty <- apply(sim[,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
  dis <- apply(sim[,1:l],MARGIN=1,FUN=pchem_d,protein2=opaa,DisMat=dismat)#distance from optimal amino acids
  dis <- apply(dis,2,sum)#sum of distances for all sites
  ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
  disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
  if(func) #plot functionality vs time
    plot(ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
  if(dist) #plot distance vs time
    plot(disfun,xlab="time",ylab="distance",main=paste("distance, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
  #return(as.numeric(tail(sim,1)[1:l])) #the sequence at the end of simulation
  return(list(sim=sim,fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
}