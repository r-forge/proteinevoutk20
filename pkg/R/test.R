gene = 3
RDatafile <- paste("~/BackupProEvo/Newton/rokas/simtree/gene",gene,".RData",sep="")
#load(RDatafile)
source("~/proteinevoutk20/pkg/R/prune.R")
obs.data <- conv(filename=fastafile,type="num") #observed data sequences
ftny.vec <- apply(obs.data,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
## find the starting and ending position of all edges in the tree
br_pos <- function(tree){
  if (is.null(attr(tree, "order")) || attr(tree, "order") !="cladewise") 
    tree <- ape:::reorder.phylo(tree,order = "cladewise")
  br <- tree$edge
  brlen <- tree$edge.length
  root <- br[1,1]
  node.pos <- vector(mode="numeric",length=max(br))
  node.pos[root] <- 0
  for(i in 1:length(brlen)){
    node.pos[br[i,2]] <- brlen[i] + node.pos[br[i,1]]
  }
  br.pos <- cbind(node.pos[br[,1]],node.pos[br[,2]])
  return(br.pos)
}
## find the path to a particular tip in a tree
## including the branches and the nodes on the path
path_to_tip <- function(tree,i){
  br <- tree$edge #all the edges
  ntips <- max(br) -  tree$Nnode
  root <- br[1,1]
  done <- FALSE
  node.path <- i #nodes on the path, from root to tip
  br.path <- NULL # indices of branches on the path to the tip
  while(!done){
    if(node.path[1]==root)
      done <- TRUE
    else{
      start <- node.path[1]
      br.ind <- which(br[,2]==start)
      br.path <- c(br.ind,br.path)
      pre <- br[br.ind,1]
      node.path <- c(pre,node.path)
    }
  }
  return(list(br.path=br.path,node.path=node.path))
}
## get functionality from simulated path, and the starting time of simulation (instead of default 0)
sim.info.t <- function(sim,opaa,t=0,s=1,beta=be,gamma=ga){
  dismat <- GM_cpv(GM_CPV,al,beta,gamma)
  l <- dim(sim)[2]-2
  steps <- dim(sim)[1]
  if(steps >= 400)
    steps <- seq(from=1,to=steps,by=steps %/% 200)
  else
    steps <- 1:steps
  
  if(length(steps) == 1)
    fty <- Ftny_protein(sim[steps,1:l],protein_op=opaa,s=s,DisMat=dismat)
  else
    fty <- apply(sim[steps,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
  time <- sim[steps,l+1]+t
  #ftyfun <- stepfun(sim[-1,l+1]+t,fty,f=0,right=FALSE)#make step functions
  return(list(fty=fty,t=time))
  #return(list(sim=sim,fty=fty,ftyfun=ftyfun))
}

plot_trace <- function(sim,model="new"){
  sim.trace <- sim$trace ## simulations on all the branches
  tree <- sim$tree ## tree on which the simulation is done
  br.pos <- br_pos(tree)
  if(model=="new")
    sim_info <- lapply(1:14, function(x) sim.info.t(sim.trace[[x]]$path,opaa[index],t=br.pos[x,1],s=s,beta,gamma))
  ## this one is for simulation under empirical model
  else
    sim_info <- lapply(1:dim(br.pos)[1], function(x) sim.info.t(sim.trace[[x]],opaa[index],t=br.pos[x,1],s=s,beta,gamma))
  ftylim <- ftyrange(sim_info)
  ntips <- length(tree$tip.label)
  par(mfrow=c(2,4))
  for(j in 1:ntips){
    tip = tree$tip.label[j]
    pathj <- path_to_tip(tree,j)
    t <- max(br.pos[pathj$br.path])
    plot(c(0,t),ftylim,type="n",bty="n",xlab="time",ylab="functionality",main=paste("tip",j),axes=FALSE,xlim=c(0,t))
    axis(1,pos=ftylim[1])
    axis(2,pos=0)
    abline(h=ftny.vec[tip])
    for(i in path_to_tip(tree,j)$br.path)
      ## plot(sim_info[[i]]$fty~sim_info[[i]]$t,pch=20,do.points=FALSE,add=TRUE)
      points(sim_info[[i]]$fty~sim_info[[i]]$t,pch=20)
  }
}