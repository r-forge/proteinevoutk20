########################################################################################################
# find the Gramma rates given the shape parameter and number of categories
discrete.gamma <- function (alpha, k) 
{
  if (k == 1) return(1)
  quants <- qgamma((1:(k - 1))/k, shape = alpha, rate = alpha)
  diff( c(0, pgamma(quants * alpha, alpha + 1),1)) * k
}
########################################################################################################
## Given s, Distance matrix, mutation rate matrix, and other parameters, find the list of 20 rate matrices corresponding to optimal aa 1:20
## same function as QAllaa, but all the matrices are scaled in this function
rate_move_mat <- function(s, DisMat, MuMat,scale.vec=rep(1,20),C=2, Phi=0.5, q=4e-7, Ne=5e6){
  fn = function(i){
    mat = mat_gen_indep(i,s,DisMat,MuMat,C,Phi,q,Ne)
    mat = mat*scale.vec[i]
  }
  mat_list = lapply(1:20,fn)
  return(mat_list)
}
rate_move <- function(protein, protein_op,s, DisMat, MuMat, scale.vec=rep(1,20),C=2, Phi=0.5, q=4e-7, Ne=5e6){
    mat = rate_move_mat(s,DisMat,MuMat,bfaa,C,Phi,q,Ne)
    vec_list = sapply(1:length(protein),function(i) {mat[[protein_op[i]]][protein[i],]},simplify="array")
    #vec_list[vec_list < 0] = 0
    return(c(vec_list))
}
##find the equilibrium frequencies for all amino acids as optimal
bfAll <- function(s,beta,gamma,GTRvec,opfreq,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  GM=GM_cpv(GM_CPV,al,beta,gamma) #Distance matrix from new weights
  mumat = aa_MuMat_form(GTRvec) #symmetric matrix for the mutation rate matrix, from the GTR matrix
  fn = function(i){
    mat = mat_gen_indep(i,s,GM,mumat,C,Phi,q,Ne)
  }
  matall = lapply(1:20,fn)

  bfall = sapply(1:20,function(i) eqmQ(matall[[i]])) # all the bf vectors, equilibrium frequencies
  bf = bfall * rep(opfreq,each=20)
  freq = apply(bf,1,sum)
  freq = freq / sum(freq)
  return(list(bfall = bfall,freq = freq))
}
########################################################################################################
#Simulation. Given the starting protein and the optimal protein
#t should be equal to expected number of substitutions
simulation <- function(protein,protein_op,t,s,DisMat,MuMat,scale.vec=rep(1,20),C=2, Phi=0.5,q=4e-7, Ne=5e6){
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site") #column names
  ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  matall = rate_move_mat(s,DisMat,MuMat,scale.vec,C,Phi,q,Ne) #Q for all optimal aa possibilities
  m = 20
  while(t_now < t){ #when current time is less than t
    vec_list = sapply(1:l,function(i) {matall[[protein_op[i]]][protein[i],]},simplify="array")
    vec_list[vec_list < 0] = 0
    rates<- c(vec_list) #moving rates of a protein to its neighbors
    lambda <- sum(rates) #total rates of moving (reaction)
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait # current time after the first transition
    ##determine where the first transition is
    if(t_now < t){
      index <- sample(1:length(rates),1,replace=T,rates) #index of the protein it moves to
      pos <- (index-1) %/% m + 1 #index of the site that changes
      rmd <- (index-1) %% m + 1 # the amino acid that site changes to
      protein[pos] <- seq(1:m)[rmd] #new protein
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  ##shift the third column up one step so that the waiting time is the time
  ## spent in the state in the first column
  path[,l+2] <- c(path[-1,l+2],NA)  
  return(list(path=path, t=t))
}
#simulation(c(1,2),c(3,4),1000,10,0.1,GM,mumat,indep=T)
#mllm(data,tree,s,beta,gamma,Q=NULL,opw=NULL,bfnu=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)

## similar to the previous function, instead of getting rate matrices form DisMat and MuMat, these are given as "matall"
simulation1 <- function(protein,protein_op,t,matall,C=2, Phi=0.5,q=4e-7, Ne=5e6){
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #column names
  ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  m = 20
  while(t_now < t){ #when current time is less than t
    vec_list = sapply(1:l,function(i) {matall[[protein_op[i]]][protein[i],]},simplify="array")
    vec_list[vec_list < 0] = 0
    rates<- c(vec_list) #moving rates of a protein to its neighbors
    lambda <- sum(rates) #total rates of moving (reaction)
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait # current time after the first transition
    ##determine where the first transition is
    if(t_now < t){
      index <- sample(1:length(rates),1,replace=T,rates) #index of the protein it moves to
      pos <- (index-1) %/% m + 1 #index of the site that changes
      rmd <- (index-1) %% m + 1 # the amino acid that site changes to
      protein[pos] <- seq(1:m)[rmd] #new protein
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  ##shift the third column up one step so that the waiting time is the time
  ## spent in the state in the first column
  path[,l+2] <- c(path[-1,l+2],NA)  
  return(list(path=path, t=t))
}

# simulation with length l, all optimal aa are the same, given by opaa.
sim_op <- function(opaa, l, t,matall,C=2, Phi=0.5,q=4e-7, Ne=5e6){
  mat = matall[[opaa]]
  bf = expm(mat*100)[1,]
  start_seq = sample(1:20,l,replace=TRUE,prob=bf)
  return(simulation1(start_seq,rep(opaa,l),t,matall,C,Phi,q,Ne))
}
##Simulation of protein sequences of length "l" on a phylogeny "tree", 
##given the ancestral sequence "rootseq", optimal amino acid sequence "protein_op",
##selection coefficient "s", here we consider s THE SAME ACROSS ALL THE SITES IN ONE GENE 
simTree <- function(tree,protein_op,s=NULL,GTRvec=NULL,alpha=al,beta=be, gamma=ga,dismat = NULL,
                    mumat=NULL,matall=NULL,scale.vec=rep(1,20),
                    rootseq=NULL,C=2,Phi=0.5,q=4e-7, Ne=5e6){
  call = match.call()
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  if (is.null(attr(tree, "order")) || attr(tree, "order") !="cladewise") 
    tree <- ape:::reorder.phylo(tree,order = "cladewise")
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree
  parent <- as.integer(edge[, 1]) #parents of the edges
  child <- as.integer(edge[, 2]) #children of the edges
  root <- as.integer(parent[!match(parent, child, 0)][1]) #root node

  tl = tree$edge.length #lengths of the edges
  sim.trace <- vector(mode="list",length=length(tl))
  m=20 #number of amino acids
  l = length(protein_op) #number of sites
  res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes
  ## get all the sub rate matrices when each amino acid is optimal
  if(is.null(matall)){
    if(is.null(dismat))
      dismat=GM_cpv(GM_CPV,alpha,beta,gamma) #Distance matrix from new weights
    if(is.null(mumat))
      mumat = aa_MuMat_form(GTRvec) #symmetric matrix for the mutation rate matrix, from the GTR matrix
    matall = rate_move_mat(s,dismat,mumat,scale.vec,C,Phi,q,Ne) #all the sub rate matrix
  }
  ## if rootseq is not given, sample it from the matrices, use equilibrium frequencies
  if(is.null(rootseq)){
    rootseq = rep(0,l) #root sequence
    bfall = lapply(1:m,function(i) eqmQ(matall[[i]])) # all the bf vectors
    for(i in 1:m){
      index = which(protein_op==i)
      rootseq[index] = sample(1:m,length(index),replace=TRUE,prob=bfall[[i]])
    }
  }
  res[root,] <- rootseq
  for(i in 1:length(tl)){
    from = parent[i] 
    to = child[i]
    ##simulation on this one branch
    sim_res = simulation1(res[from,],protein_op,tl[i],matall,C,Phi,q,Ne)
    res[to,] = as.numeric(tail(sim_res$path,1)[1:l])
    sim.trace[[i]] <- sim_res
  }
  res = AA[res]
  res = matrix(res,ncol=l)
  k = length(tree$tip)
  label = c(tree$tip, as.character((k+1):nNodes))
  rownames(res)=label 
  res = res[ tree$tip, , drop=FALSE]
  return(list(data=res,trace=sim.trace,tree=tree,call=call))
}
#simTree(tree,5,protein_op=Protein_op[1:5],m=20,s=0.3,Nu_vec,al,be,ga,q=4e-7,Ne=5e6,Root[1:5])
simTreeEmp <- function(tree,l=100,rootseq=NULL,Q=NULL,bf=NULL,inv=0,rate=1,k=1,model="USER"){
  call = match.call()
  if(!is.null(model) && model!="USER"){
    model <- match.arg(model,get(".aamodels",environment(pml)))
    getModelAA(model,bf=is.null(bf),Q=is.null(Q))
  }
  m=20 #number of amino acids
  if(is.null(bf)) bf = rep(1/m,m)
  if(is.null(Q)) Q = rep(1,m*(m-1)/2)
  
  if(!is.binary.tree(tree)) stop("error: the input phylogeny is not binary tree!")
  if (is.null(attr(tree, "order")) || attr(tree, "order") !="cladewise") 
    tree <- ape:::reorder.phylo(tree,order = "cladewise")
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree
  parent <- as.integer(edge[, 1]) #parents of the edges
  child <- as.integer(edge[, 2]) #children of the edges
  root <- as.integer(parent[!match(parent, child, 0)][1]) #root node
  if(is.null(rootseq)) rootseq=sample(1:20,l,replace=TRUE,prob=bf)
  
  tl = tree$edge.length #lengths of the edges
  sim.trace <- vector(mode="list",length=length(tl))

  l = length(rootseq) #number of sites
  res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes
  res[root,] <- rootseq
  for(i in 1:length(tl)){
    from = parent[i] 
    to = child[i]
    ##simulation on this one branch
    sim_res = simulationQ(protein=res[from,],t=tl[i],Q=Q,bf=bf,inv=inv,rate=rate,k=k)
    res[to,] = as.numeric(tail(sim_res$path,1)[1:l])
    sim.trace[[i]] <- sim_res
  }
  res = AA[res]
  res = matrix(res,ncol=l)
  k = length(tree$tip)
  label = c(tree$tip, as.character((k+1):nNodes))
  rownames(res)=label 
  res = res[ tree$tip, , drop=FALSE]
  return(list(data=res,trace=sim.trace,tree=tree,call=call))
}
#######################################################################
## simulation along the tree, given data, best model and tree
sim_emp <- function(data,best_emp_model){
  nsites <- length(attr(data,"index"))
  model <- best_emp_model$model
  is.G <- "G" %in% model #gamma?
  is.F <- "F" %in% model #empirical frequencies?
  is.I <- "I" %in% model #invariant sites included?
  bf = NULL
  k = 1
  shape = 1
  if(is.F) bf=findBf2(data)
  if(is.I) inv=best_emp_model$inv
  if(is.G) {shape = best_emp_model$shape 
            k=4}
  sim.emp <- simTreeEmp(best_emp_model$tree,l=nsites,bf=bf,rate=shape,k=k,model=model[1])
  return(sim.emp)
}
##################################################################################################
## simulation under the empirical models
##################################################################################################
## simulation given rate matrix (or lower triangular part) Q, proportion of invariant sites inv,
## Gamma rate, and number of categories k, base frequencies bf.
simulationQ <- function(protein,t,Q=NULL,bf=NULL,inv=0,rate=1,k=1){
  if(is.null(bf)) bf <- rep(1/20,20)
  if(is.null(Q)) Q <- rep(1,190)
  if(is.matrix(Q)) Qmat <- Q
  else if(is.vector(Q))  
    Qmat <- mat_form_lowtriQ(Q,bf,byrow=TRUE) #transition rate matrix
  l <- length(protein) #number of sites
  
  inv_runif <- runif(l) #sample the sites that are invariant
  inv_site <- which(inv_runif < inv) #this could have length 0 if inv is very small
  grates <- discrete.gamma(rate,k) #rates in different categories for gamma distribution
  var_site <- setdiff(1:l,inv_site) #indices of variant sites
  protein_inv <- protein[inv_site] #invariant sites
  protein_var <- protein[var_site] #simulation on the variant sites only, then attach the invariant sites to the variant sites
  lvar <- length(var_site) #number of variant sites
  var_cat <- sample(1:k,lvar,replace=T,prob=rep(1,k)) #assign gamma categories for variant sites
  
  sprotein <- protein_var #starting protein
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #column names
  ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  m = 20
  rwait <- -diag(Qmat) #exponential rate of waiting, diagonal of Qmat
  while(t_now < t){ #when current time is less than t
    rates <- Qmat[sprotein,] * grates[var_cat]
    rates <- c(t(rates)) #moving rates of a protein to its neighbors
    rates[rates < 0] <- 0    
    #rates <- rwait[sprotein]*grates[var_cat] #take into account of gamma rates
    lambda <- sum(rates) #total rates of moving (reaction)
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait # current time after the first transition
    ##determine where the first transition is
    if(t_now < t){
      index <- sample(1:length(rates),1,replace=T,rates) #index of the protein it moves to
      pos <- (index-1) %/% m + 1 #index of the site that changes
      rmd <- (index-1) %% m + 1 # the amino acid that site changes to
      sprotein[pos] <- seq(1:m)[rmd] #new protein
      protein[var_site] <- sprotein
      protein[inv_site] <- protein_inv
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  path[,l+2] <- c(path[-1,l+2],NA)  
  return(list(path=path,t=t))
}
## simulation according to empirical models,l is the length of simulated sequence, rootseq is the starting sequence
## if none of Q and bf is specified, then they are got from the model
## this is just along one branch with given length, not along a tree
simAA <- function(l=100,rootseq=NULL,t=1,Q=NULL,bf=NULL,inv=0,rate=1,k=1,model="USER"){
  if(!is.null(model)){
    model <- match.arg(model,get(".aamodels",environment(pml)))
    if(model!="USER") getModelAA(model,bf=is.null(bf),Q=is.null(Q))
  }
  lbf = 20
  if(is.null(bf)) bf = rep(1/lbf,lbf)
  if(is.null(Q))  Q = rep(1,lbf*(lbf-1)/2)
  if(is.null(rootseq)) rootseq=sample(1:20,l,replace=TRUE,prob=bf)
  return(simulationQ(protein=rootseq,t=t,Q=Q,bf=bf,inv=inv,rate=rate,k=k))
}
##################################################################################################
## likelihood of tree, given observed sequences, own code, for check
##################################################################################################
## find the site likelihood given tree, data, and the lower triangular part of Q
## only works for data of one site, with integers as states
# ll_site_lowQ <- function(tree,data,Q,bf=rep(1/20,20),g=1){
#   ##If the given tree is not rooted and binary, then throw error and exit
#   if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
#   m = 20
#   ##if the base frequencies are not specified, do a uniform distribution
#   if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
#   Qmat <- mat_form_lowtriQ(Q,bf,byrow=TRUE)
#   
#   tree <- ape:::reorder.phylo(tree,"p") #reorder the tree in pruningwise order
#   edge = tree$edge #edges
#   nNodes = max(edge) #number of nodes in the tree (including tips)
#   probvec = matrix(NA,nNodes,m) #probability of gettting different states at nodes that evolve to the current sequences
#   parent <- as.integer(edge[, 1]) #parents of the edges
#   child <- as.integer(edge[, 2]) #children of the edges
#   root <- as.integer(parent[!match(parent, child, 0)][1])  
#   tip <- as.integer(child[!match(child,parent,0)])
#   
#   init.tip <- function(x){ #initiate the vector for the tips, 1 for the tip state, 0 otherwise
#     vec <- rep(0,m)
#     vec[data[x]] <- 1
#     vec
#   }
#   probvec[1:length(tip),] <- t(sapply(1:length(tip),init.tip)) #all tips
#   tl = tree$edge.length #lengths of the edges
#   P <- getPm(tl,Qmat,1)
#   for(i in 1:tree$Nnode){ #for each interior node calculate the probability vector of observing 1 of 20 states
#     from = parent[2*i] #parents
#     to = child[(2*i-1):(2*i)] #direct descendents
#     v.left <- P[[1,2*i-1]]
#     v.right <- P[[1,2*i]]
#     probvec[from,] <- as.vector((v.left%*%probvec[to[1],])*(v.right%*%probvec[to[2],])) #pruning, vector form
#     check.sum <- sum(probvec[from,])
#     if(check.sum==0) #probability is very very low
#       warning("numerical overflow",immediate.=TRUE)
#   }
#   return(probvec[root,])
# }
##################################################################################################
## given the simulation sequences with times of transitions, find the functionalities of each protein 
## on the path, and the distance of them from the optimal aa, as well as the stepfunction formed by them
## The simulated sequence data will have at least 2 sequences, that's when no site encounters
## any substitution. t: starting time of simulation
sim.info <- function(sim,opaa,obsaa=NULL,ratio=TRUE,t=0,s=1,beta=be,gamma=ga,max.step=200,fty=TRUE,dist=TRUE){
  dismat <- GM_cpv(GM_CPV,al,beta,gamma)
  l <- dim(sim)[2]-2 # number of sites
  ### if there are too many steps, choose less than max.step(default 200) to shorten the running time
  steps <- dim(sim)[1]
  if(steps >= 400)
    steps <- seq(from=1,to=steps,by=steps %/% max.step)
  else
    steps <- 1:steps
  if(fty){
    if(length(steps) == 1)
      fty <- Ftny_protein(sim[steps,1:l],protein_op=opaa,s=s,DisMat=dismat)
    else
      fty <- apply(sim[steps,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
  }
  time <- sim[steps,l+1]+t #time of each step
  #ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
  if((!is.null(obsaa))&&dist){
    if(!ratio){
      dis <- apply(sim[steps,1:l],MARGIN=1,FUN=pchem_d,protein2=obsaa,DisMat=dismat)#distance from optimal amino acids
      dis <- -apply(dis,2,mean)#average distance for all sites, opposite sign
      #disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
      # return(list(fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
      return(list(fty=fty,dis=dis,t=time))
    }
    else{
      dis <- sapply(steps, function(x) sum(sim[x,1:l]==obsaa)/l)
      return(list(fty=fty,dis=dis,t=time))
    }
  }
    
  else
    return(list(fty=fty,t=time))
}

## find the range of the functionalities from a list of objects from sim.info
## uses $fty component
ftyrange <- function(sim_info){
  nsim <- length(sim_info)
  res <- NULL
  res <- cbind(sapply(1:nsim, function(x) range(sim_info[[x]]$fty)))
  res <- c(min(res[1,]),max(res[2,]))
  return(res)
}
## find the range of distance between observed sequence form a list of objs from sim.info
## uses $dis component
disrange <- function(sim_info){
  nsim <- length(sim_info)
  res <- NULL
  res <- cbind(sapply(1:nsim, function(x) range(sim_info[[x]]$dis)))
  res <- c(min(res[1,]),max(res[2,]))
  return(res)
}
  ## given amino acid model, assign sub matrix Q and base frequencies bf to parent frame
getModelAA <- function(model, bf=TRUE, Q=TRUE){
  model <- match.arg(eval(model),get(".aamodels",environment(pml)))
  tmp = get(paste(".", model, sep=""),environment(pml))
  if(Q) assign("Q", tmp$Q, envir=parent.frame())
  if(bf) assign("bf", tmp$bf, envir=parent.frame())
}

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
get_trace <- function(sim,s,beta,gamma,ftny.vec,ratio=TRUE){
  sim.trace <- sim$trace ## simulations on all the branches
  l <- length(sim.trace) ## number of traces (branches)
  tree <- sim$tree ## tree on which the simulation is done
  br.pos <- br_pos(tree) ##edges, with entries as the position of nodes
  ntips <- length(tree$tip.label) #number of tips
  
  trace.info <- vector("list",length=ntips)
  for(j in 1:ntips){
    tip = tree$tip.label[j] #name of the tip
    obs.ftny <- ftny.vec[tip]
    pathj <- path_to_tip(tree,j)
    pathj.bgn <- br.pos[pathj$br.path,1] #begining position of the branches on the path 
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
    trace.info[[j]]$tip <- tip
  }
  return(trace.info)
}
plot_trace <- function(trace.info,obs.ftny,zoom=FALSE,plotftny=TRUE,plotdis=FALSE,label="gene",gene=1){
  if(plotftny){
    par(mfrow=c(2,4))
    ftylim <- ftyrange(trace.info)
    if(!zoom)
      ftylim[2] <- 1
    ftylim<- range(c(ftylim,obs.ftny))
    for(i in 1:length(trace.info)){
      sim_info <- trace.info[[i]]
      t <- max(sim_info$t)
      plot(c(0,t),ftylim,type="n",bty="n",xlab="time",ylab="functionality",
           main=paste(label, gene, ",", sim_info$tip),axes=FALSE,xlim=c(0,t))
      axis(1,pos=ftylim[1])
      axis(2,pos=0)
      abline(h=obs.ftny[sim_info$tip],col="blue")
      points(sim_info$fty~sim_info$t,pch=20)
    }
  }
  ## plot distance
  if(plotdis){
    par(mfrow=c(2,4))
    dislim <- disrange(trace.info)
    dislim[2] <- 1
    for(i in 1:length(trace.info)){
      sim_info <- trace.info[[i]]
      t <- max(sim_info$t)
      plot(c(0,t),dislim,type="n",bty="n",xlab="time",ylab="similarity",
           main=paste(label, gene, ",", sim_info$tip),axes=FALSE,xlim=c(0,t))
      axis(1,pos=dislim[1])
      axis(2,pos=0)
      points(sim_info$dis~sim_info$t,pch=20)
    }
  }
}