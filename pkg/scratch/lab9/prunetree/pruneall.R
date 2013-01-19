system("sed -e '11d;12d' ~/proteinevoutk20/pkg/Data/gene1.fasta > ~/proteinevoutk20/pkg/scratch/lab9/prunetree/gene1p.fasta")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
## import data on 7 tips and 8 tips
datap = conv("~/proteinevoutk20/pkg/scratch/lab9/prunetree/gene1p.fasta","phyDat")
data = conv("~/proteinevoutk20/pkg/Data/gene1.fasta","phyDat")
dataPnum <- matrix(unlist(datap),nrow=7,byrow=TRUE) #7-tip data in numbers stored in matrix
datanum <- conv("~/proteinevoutk20/pkg/Data/gene1.fasta")
datanum <- matrix(unlist(datanum),nrow=8,byrow=TRUE) #8-tip data in numbers stored in matrix
## starting tree with Scas pruned
tre = read.tree("~/proteinevoutk20/pkg/Data/prune.tre")
res = mllm1(datap,tre,0.1,be,ga,Q=NU_VEC)
res$ll$loglik
## optimal parameters with 7 tip data
res_op <- optim.mllm(res,optQ=T,optBranch=T,optsWeight=T,optOpw=FALSE,
                     control=list(epsilon=1e-08,hmaxit=500,htrace=1,print_level=0,maxeval="50"))
res_op$ll$loglik
s = res_op$s
beta=res_op$GMweights[2]
gamma = res_op$GMweights[3]
Q = res_op$Q
treep = res_op$tree
index = attr(datap,"index")
opaa = res_op$ll$opaa[index]
opaap = res_op$ll$opaa
bfaa = res_op$bfaa
dismat = res_op$dismat
mumat = res_op$mumat
tree = ape:::reorder.phylo(ROKAS_TREE,"p")
## optimize branch lengths on 8-tip tree with parameters just found, and 8-tip data
br <- optim.br(data,tree,maxeval="30",print_level=1,s=s,beta=beta,gamma=gamma,Q=Q,opaa=opaa)
br
tree$edge.length <- br$solution
tree$edge.length[10] <- 0
brlen <- br$solution[10]

nr <- attr(datap,"nr") #number of different patterns in 7-tip data
# for site i, find the probability (likelihood) of all 20 states at that site
sitei <- function(i){
  datai <- dataPnum[,i] #ith distinct pattern in the truncated(pruned) data                                                                                     
  addj <- function(j){
    datai_j <- append(datai,values=j,after=5) #Add Scas as a tip but with external branch length equal to 0                                                     
    return(ll_site(tree,datai_j,optimal=opaap[i],s=s,Q=Q,alpha=al,beta=beta,gamma=gamma,bf=bfaa)) #loglikelihood at this site with Scas has state j             
  }
  sitelike <- sapply(1:20,addj) # all 20 possible states at site i                                                                                              
  return(sitelike)
}

siteprob <- lapply(1:nr,sitei) # do all sites
siteprob <- simplify2array(siteprob,higher=FALSE) #convert result from list to matrix, 20 by "nr" matrix
siteprop <- apply(siteprob,2,function(x) x/sum(x)) #normalize each site prob so that the sum is 1, proportion to be used in simulation 
siteprop <- siteprop[,index]

nsim <- 100
##simulation under our new model
sim <- vector(mode="list",length=nsim)
simWag <- vector(mode="list",length=nsim)
roots <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  root <- sapply(1:length(index),function(x) sample(20,1,replace=TRUE,prob=siteprop[,x]))
  roots[[i]] <- root
  sim[[i]] <- sim.New(s=s,t=brlen,root=root,opaa=opaa,beta=beta,gamma=gamma,bfaa=bfaa)
  simWag[[i]] <- sim.Wag(t=brlen,protein=root,bf=bfaa,opaa=opaa,s=s,dismat=dismat)
}
##################################################################################################
## analysis of the simulated data under both models
##################################################################################################
pdf("~/proteinevoutk20/pkg/Plot/prunetree/gene1.pdf")
##functionalities based on the parameters estimated earlier
start.ftny.vec <-  sapply(1:nsim,function(i) head(sim[[i]]$fty,1))
end.ftny.vec <- sapply(1:nsim,function(i) tail(sim[[i]]$fty,1))
plot.density(density(end.ftny.vec),xlab="Functionality",ylab="density",main="Our model")
abline(v=Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat))
##Grantham distances from optimal aa's
start.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=roots[[i]],protein2=opaa,DisMat=GM)))
end.dis.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(sim[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
plot.density(density(end.dis.vec),xlab="Avg Distance from optima",ylab="density",main="Our model")
abline(v=mean(pchem_d(datanum[6,],opaa,DisMat=dismat)))
##Grantham distances from observed sequence
end.dis.obs.vec <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(sim[[x]]$sim[,1:length(index)],1)),protein2=datanum[6,],DisMat=GM)))
plot.density(density(end.dis.obs.vec),xlab="Avg Distance from observed seq",ylab="density",main="Our model")
##################################################################################################
##WAG model
end.ftny.vec.wag <- sapply(1:nsim,function(i) tail(simWag[[i]]$fty,1))
plot.density(density(end.ftny.vec.wag),xlab="Functionality",ylab="density",main="WAG model")
abline(v=Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat))
##Grantham distances from optimal aa's
end.dis.vec.wag <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(simWag[[x]]$sim[,1:length(index)],1)),protein2=opaa,DisMat=GM)))
plot.density(density(end.dis.vec.wag),xlab="Avg Distance from optima",ylab="density",main="WAG model")
abline(v=mean(pchem_d(datanum[6,],opaa,DisMat=dismat)))
##Grantham distances from observed sequence
end.dis.obs.vec.wag <- sapply(1:nsim,function(x) mean(pchem_d(protein1=as.numeric(tail(simWag[[x]]$sim[,1:length(index)],1)),protein2=datanum[6,],DisMat=GM)))
plot.density(density(end.dis.obs.vec.wag),xlab="Avg Distance from observed seq",ylab="density",main="WAG model")
##################################################################################################
dev.off()