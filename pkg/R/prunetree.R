load("~/proteinevoutk20/pkg/scratch/lab9/gene7p/gene7p.RData")
source("~/proteinevoutk20/pkg/R/main.R")
## the objects for the pruned data/tree has suffix "p" indicating that it's for the pruned
## use the paraemters found with 7 tips, find the optimal branch lengths with 8 tips
data = res_op$data
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

tree <- ROKAS_TREE
#br <- optim.br(gene7,tree,maxeval="3000",print_level=1,s=s,beta=beta,gamma=gamma,Q=Q,opaa=opaa)

#load("~/proteinevoutk20/pkg/scratch/lab9/gene7p/gene7p.RData")
tree <- ape:::reorder.phylo(tree,"p")
tree$edge.length[10] <- 0 # set the edge leading to Scas to 0
nr <- attr(gene7p,"nr")
opaap <- res_op$ll$opaa # optimal amino acids at all different site patterns, instead of all sites
Qall <- res_op$Qall
#for each site pattern:

sitei <- function(i){
  datai <- subset(gene7p,1:7,i,site.pattern=FALSE) #ith distinct pattern in the truncated(pruned) data
  addj <- function(j){
    datai$Scas <- j #Add Scas as a tip but with external branch length equal to 0
    return(mllm1(datai,tree,Qall=Qall,bfaa=bfaa,opaa=opaap[i])$ll$loglik) #loglikelihood at this site with Scas has state j
  }
  sitelike <- sapply(1:20,addj) # all 20 possible states at site i
  return(sitelike)
}

siteprob <- mclapply(1:nr,sitei)
siteprob <- simplify2array(siteprob,higher=FALSE)
siteprob <- exp(siteprob) #probability of observing state i at site j
siteprop <- apply(siteprob,2,function(x) x/sum(x)) #normalize each site prob so that the sum is 1, proportion to be used in simulation
save(siteprop,siteprob,file="proteinevoutk20/pkg/scratch/lab9/gene7p/siteprob.RData",compress=TRUE)

siteprop <- siteprop[,index]
root<- sapply(1:length(index),function(x) sample(20,1,replace=TRUE,prob=siteprop[,x]))
# in Amino Acid: 
# sample <- AA[sample]
brlen <- br$solution[10]
nsim <- 100
sim <- vector(mode="list",length=nsim)
sim_start <- plot.sim(s=s,t=brlen,root=root,opaa=opaa,beta=beta,gamma=gamma,bfaa=bfaa,dist=FALSE,add=F)
roots <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  root<- sapply(1:length(index),function(x) sample(20,1,replace=TRUE,prob=siteprop[,x]))
  roots[[i]] <- root
  sim[[i]] <- plot.sim(s=s,t=brlen,root=root,opaa=opaa,beta=beta,gamma=gamma,bfaa=bfaa,dist=FALSE,add=T)
}
ftny.vec <- sapply(1:nsim,function(i) tail(sim[[i]]$fty,1))
plot.density(density(ftny.vec))
abline(v=Ftny_protein(protein=gene7$Scas[index],protein_op=opaa,s=s,DisMat=dismat))

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