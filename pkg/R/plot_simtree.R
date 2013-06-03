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

plot_trace <- function(sim,plottip=TRUE){
  sim.trace <- sim$trace ## simulations on all the branches
  l <- length(sim.trace) ## number of traces (branches)
  tree <- sim$tree ## tree on which the simulation is done
  br.pos <- br_pos(tree) ##edges, with entries as the position of nodes
  sim_info <- lapply(1:l, function(x) sim.info(sim.trace[[x]]$path,opaa=opaa[index],t=br.pos[x,1],s=s,beta=beta,gamma=gamma))
  
  ftylim <- ftyrange(sim_info) #range of the functionalities
  ntips <- length(tree$tip.label) #number of tips
  par(mfrow=c(2,4))
  for(j in 1:ntips){
    tip = tree$tip.label[j] #name of the tip
    obs.ftny <- ftny.vec[tip]
    pathj <- path_to_tip(tree,j)
    t <- max(br.pos[pathj$br.path])
    if(plottip)
      ftylim.tip <- range(c(ftylim,obs.ftny))
    else 
      ftylim.tip <- ftylim
    plot(c(0,t),ftylim.tip,type="n",bty="n",xlab="time",ylab="functionality",main=paste("gene", gene, ",", tip),axes=FALSE,xlim=c(0,t))
    axis(1,pos=ftylim.tip[1])
    axis(2,pos=0)
    abline(h=obs.ftny,col="blue")
    for(i in path_to_tip(tree,j)$br.path)
      ## plot(sim_info[[i]]$fty~sim_info[[i]]$t,pch=20,do.points=FALSE,add=TRUE)
      points(sim_info[[i]]$fty~sim_info[[i]]$t,pch=20)
  }
}