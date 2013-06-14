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
  ntips <- length(tree$tip.label)
  parents <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parents[!match(parents, child, 0)][1]) 
  
  done <- FALSE
  node.path <- i #nodes on the path, from root to tip
  br.path <- NULL # indices of branches on the path to the tip
  while(!done){
    if(node.path[1]==root)
      done <- TRUE
    else{
      start <- node.path[1] #start tracing, go backwards
      br.ind <- which(br[,2]==start) #find the index of branch that ends in "start"
      br.path <- c(br.ind,br.path)
      pre <- br[br.ind,1] #next node on the path when going backwards
      node.path <- c(pre,node.path)
    }
  }
  return(list(br.path=br.path,node.path=node.path))
}

plot_trace <- function(sim,ratio=TRUE,plotftny=TRUE,plotdis=FALSE,plottip=TRUE){
  sim.trace <- sim$trace ## simulations on all the branches
  l <- length(sim.trace) ## number of traces (branches)
  tree <- sim$tree ## tree on which the simulation is done
  br.pos <- br_pos(tree) ##edges, with entries as the position of nodes
  #sim_info <- lapply(1:l, function(x) sim.info(sim.trace[[x]]$path,opaa=opaa[index],t=br.pos[x,1],s=s,beta=beta,gamma=gamma))
  
  #ftylim <- ftyrange(sim_info) #range of the functionalities
  ntips <- length(tree$tip.label) #number of tips
  
  trace.info <- vector("list",length=ntips)
  for(j in 1:ntips){
    tip = tree$tip.label[j] #name of the tip
    obs.ftny <- ftny.vec[tip]
    pathj <- path_to_tip(tree,j)
    pathj.bgn <- br.pos[pathj$br.path,1]
    ########################################
    ## paste traces together for a certain tip
    tracej <- sim.trace[[pathj$br.path[1]]]$path
    if(length(pathj$br.path) > 1){
      for(k in 2:length(pathj.bgn)){
        trace.add <- sim.trace[[pathj$br.path[[k]]]]$path
        trace.add[,"Time Now"] <- trace.add[,"Time Now"] + pathj.bgn[k]
        tracej <- rbind(tracej,trace.add[-1,])
      }
    }
    ########################################
    t <- max(br.pos[pathj$br.path,])
    sim_info <- sim.info(sim=tracej,opaa=opaa[index],obsaa=datanum[tip,],ratio=ratio,t=0,s=s,beta=beta,gamma=gamma,fty=T,dist=T)
    trace.info[[j]] <- sim_info
    if(plotftny){
      par(mfrow=c(2,4))
      ftylim <- range(sim_info$fty)
      ftylim[2] <- 1
      ftylim<- range(c(ftylim,obs.ftny))
      
      plot(c(0,t),ftylim,type="n",bty="n",xlab="time",ylab="functionality",main=paste("gene", gene, ",", tip),axes=FALSE,xlim=c(0,t))
      axis(1,pos=ftylim[1])
      axis(2,pos=0)
      abline(h=obs.ftny,col="blue")
      points(sim_info$fty~sim_info$t,pch=20)
    }
    ## plot distance
    if(plotdis){
      par(mfrow=c(2,4))
      dislim <- range(sim_info$dis)
      dislim[2] <- 1
      plot(c(0,t),dislim,type="n",bty="n",xlab="time",ylab="similarity",main=paste("gene", gene, ",", tip),axes=FALSE,xlim=c(0,t))
      axis(1,pos=dislim[1])
      axis(2,pos=0)
      points(sim_info$dis~sim_info$t,pch=20)
    }
  }
  return(trace.info)
}