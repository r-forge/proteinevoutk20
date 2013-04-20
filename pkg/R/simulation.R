## source the R file that contains required functions
#source("~/proteinevoutk20/pkg/R/pphyproevo.R")
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
rate_move_mat <- function(s, DisMat, MuMat,bfaa = rep(1/20,20),C=2, Phi=0.5, q=4e-7, Ne=5e6){
  fn = function(i){
    mat = mat_gen_indep(i,s,DisMat,MuMat,C,Phi,q,Ne)
    mat = scaleQ(mat,bfaa)
  }
  mat_list = lapply(1:20,fn)
  return(mat_list)
}
rate_move <- function(protein, protein_op,s, DisMat, MuMat, bfaa = rep(1/20,20),C=2, Phi=0.5, q=4e-7, Ne=5e6){
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

  bfall = sapply(1:20,function(i) expm(matall[[i]]*100)[1,]) # all the bf vectors
  bf = bfall * rep(opfreq,each=20)
  freq = apply(bf,1,sum)
  freq = freq / sum(freq)
  return(list(bfall = bfall,freq = freq))
}
########################################################################################################
#Simulation. Given the starting protein and the optimal protein
#t should be equal to expected number of substitutions
simulation <- function(protein,protein_op,t,s,DisMat,MuMat,bfaa=rep(1/20,20),C=2, Phi=0.5,q=4e-7, Ne=5e6){
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site") #column names
  ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  matall = rate_move_mat(s,DisMat,MuMat,bfaa,C,Phi,q,Ne) #Q for all optimal aa possibilities
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
  path
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
  return(list(path=path,start_seq=protein,op_seq=protein_op, t=t))
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
simTree <- function(tree,protein_op,s,GTRvec,alpha=al,beta=be, gamma=ga,mumat=NULL,bfaa=rep(1/20,20),
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
  m=20
  l = length(protein_op) #number of sites
  res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes

  GM=GM_cpv(GM_CPV,alpha,beta,gamma) #Distance matrix from new weights
  if(is.null(mumat)){
    mumat = aa_MuMat_form(GTRvec) #symmetric matrix for the mutation rate matrix, from the GTR matrix
    #mumat = sym.to.Q(mumat,bfaa) #from symmetric matrix to rate matrix
  }
  matall = rate_move_mat(s,GM,mumat,bfaa,C,Phi,q,Ne) #all the sub rate matrix
  
  if(is.null(rootseq)){
    rootseq = rep(0,l) #root sequence
    bfall = lapply(1:m,function(i) expm(matall[[i]]*100)[1,]) # all the bf vectors
    ## figure out the starting sequence
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
  }
  res = AA[res]
  res = matrix(res,ncol=l)
  k = length(tree$tip)
  label = c(tree$tip, as.character((k+1):nNodes))
  rownames(res)=label 
  res = res[ tree$tip, , drop=FALSE]
  return(list(data=res,tree=tree,optimal=protein_op,rootseq=rootseq,s=s,Q=GTRvec,GMweights=c(alpha,beta,gamma),bfaa=bfaa,call=call))
}
#simTree(tree,5,protein_op=Protein_op[1:5],m=20,s=0.3,Nu_vec,al,be,ga,q=4e-7,Ne=5e6,Root[1:5])

##################################################################################################
## simulation under the new model and WAG model
##################################################################################################
## Read the WAG matrix (lower triangular part)
WagMat <- scan("~/proteinevoutk20/pkg/Data/AAmodel/wag.txt",nlines=19)
#Qwag <- mat_form_lowtriQ(Q=WagMat,bf=bfaa)
ll_site_lowQ <- function(tree,data,Q,bf=rep(1/20,20),g=1){
  ##If the given tree is not rooted and binary, then throw error and exit
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  m = 20
  ##if the base frequencies are not specified, do a uniform distribution
  if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
  Qmat <- mat_form_lowtriQ(Q,bf)
  
  tree <- ape:::reorder.phylo(tree,"p") #reorder the tree in pruningwise order
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree (including tips)
  probvec = matrix(NA,nNodes,m) #probability of gettting different states at nodes that evolve to the current sequences
  parent <- as.integer(edge[, 1]) #parents of the edges
  child <- as.integer(edge[, 2]) #children of the edges
  root <- as.integer(parent[!match(parent, child, 0)][1])  
  tip <- as.integer(child[!match(child,parent,0)])
  
  init.tip <- function(x){ #initiate the vector for the tips, 1 for the tip state, 0 otherwise
    vec <- rep(0,m)
    vec[data[x]] <- 1
    vec
  }
  probvec[1:length(tip),] <- t(sapply(1:length(tip),init.tip)) #all tips
  tl = tree$edge.length #lengths of the edges
  P <- getPm(tl,Qmat,1)
  for(i in 1:tree$Nnode){ #for each interior node calculate the probability vector of observing 1 of 20 states
    from = parent[2*i] #parents
    to = child[(2*i-1):(2*i)] #direct descendents
    v.left <- P[[1,2*i-1]]
    v.right <- P[[1,2*i]]
    probvec[from,] <- as.vector((v.left%*%probvec[to[1],])*(v.right%*%probvec[to[2],])) #pruning, vector form
    check.sum <- sum(probvec[from,])
    if(check.sum==0) #probability is very very low
      warning("numerical overflow",immediate.=TRUE)
  }
  return(probvec[root,])
  #return(list(ll=max(probvec[root,]),root=which.max(probvec[root,]))) #with the corresponding root returned 
  #return(max(probvec[root,])) #just the value
}
## simulation given rate matrix, 
simulationQ <- function(protein,t,Q,bf=rep(1/20,20)){
  Qmat <- mat_form_lowtriQ(Q,bf)
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #column names
  ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  m = 20
  while(t_now < t){ #when current time is less than t
    rates <- Qmat[protein,]
    rates <- c(t(rates)) #moving rates of a protein to its neighbors
    rates[rates < 0] <- 0
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
  return(path)
}
##################################################################################################
##result includes the functionality of the protein at each step, mean distance from the optimal protein
## and the corresponding step functions of fty and distance vectors
## simulation based on WAG model
sim.Wag <- function(t=1,protein,bf,opaa,s,dismat){
  sim <- simulationQ(protein=protein,t=t,Q=WagMat,bf=bf)
  l <- dim(sim)[2] - 2
  fty <- apply(sim[,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
  dis <- apply(sim[,1:l],MARGIN=1,FUN=pchem_d,protein2=opaa,DisMat=dismat)#distance from optimal amino acids
  dis <- apply(dis,2,mean)#average distance for all sites
  ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
  disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
  #plot(ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
  return(list(sim=sim,fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
}
## simulation based on SAC model, given s, opaa, beta, gamma, bfaa,and starting sequence on a single branch with length t.
sim.New <- function(s=1,t=10,root,opaa,beta,gamma,bfaa){
  dismat = GM_cpv(GM_CPV,al,beta,gamma)
  mumat = aa_MuMat_form(res_op$Q)
  sim <- simulation(root,opaa,t=t,s=s,DisMat=dismat,MuMat=mumat,bfaa=bfaa) #simulation
  l <- dim(sim)[2] - 2
  fty <- apply(sim[,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
  dis <- apply(sim[,1:l],MARGIN=1,FUN=pchem_d,protein2=opaa,DisMat=dismat)#distance from optimal amino acids
  dis <- apply(dis,2,mean)#average distance for all sites
  ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
  disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
  #if(func) #plot functionality vs time
  #plot(ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
  #if(dist) #plot distance vs time
  #plot(disfun,xlab="time",ylab="distance",main=paste("distance, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
  #return(as.numeric(tail(sim,1)[1:l])) #the sequence at the end of simulation
  return(list(sim=sim,fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
}

##################################################################################################
## wrapper for seq-gen (seq-gen already installed on the computer)
## opts is the string of options one would use in the command "seq-gen"
## opts: [-m MODEL] [-l LENGTH] [-n NUM OF DATASETS] [-p NUM OF PARTITIONS] [-s BRANCH LENGTH SCALING] [-d TOTAL TREE SCALE]
## [-k ANCESTRAL SEQUENCE] [-c RATES FOR CODON POSITION HETEROGENEITY, 3 numbers] [-a SHAPE OF GAMMA DISN] [-g GAMMA CATEGORIES] 
## [-i PROPORTION OF INVARIANT SITES] [-t TRANSITION TRANSVERSION RATIO] [-r general rate matrix, 6 numbers][-f e BASE FREQUENCIES]
## [-z SEED FOR RANDOM NUMBER GENERATOR] [-o prn OUTPUT FILE FORMAT: phylip/relaxed phylip/nexus]
## [-w ar WRITE ADDITIONAL INFO: ancestral seq/rate for each site] [-q QUIET]
##################################################################################################
seq_gen <- function(opts, inputfile, outputfile){
  command <- paste("seq-gen", opts,"<",inputfile, ">",outputfile,sep=" ")
  system(command)
}
## opts <- "-mWAG -l144 -a3.514 -i0.193 -f0.03559028,0.05121528,0.02777778,0.01649306,0.01388889,0.03472222,0.05121528,0.06944444,0.02777778,0.09375000,0.14409722,0.05034722,0.01909722,0.05815972,0.02864583,0.07552083,0.04687500,0.02343750,0.04600694,0.08593750 -on -war"
## seq_gen(opts,"gene7.tree","gene7sim.txt")