## source the R file that contains required functions
source("~/proteinevoutk20/pkg/R/pphyproevo.R")


## Given a protein, the optimal protein, and the selection coefficient, mutation matrix and distance matrix,
## Return a matrix with each column as the rates of moving from current state to other states for the corresponding site
## Change: concatenate all the columns in the matrix to a vector, and return as result
## possible optimization: if two sites have the same current state and optimal state, the computation is exactly the same, only need once

rate_move <- function(protein, protein_op, s, DisMat, MuMat, m=20, C=2, Phi=0.5, q=4e-7, Ne=1.36e7){
  mat_move <- sapply(1:length(protein), function(i) {mat_gen_indep(protein_op[i],s,DisMat, MuMat, m, C, Phi, q, Ne)[protein[i],]}, simplify="array")
  mat_move[mat_move < 0] <- 0
  return(c(mat_move))
}

#Simulation. Given the starting protein and the optimal protein
#t: running time of the chain
#MuMat: matrix of mutation rates between amino acids
simulation <- function(protein,protein_op,t,m,s,DisMat,MuMat, C=2, Phi=0.5,q=4e-7, Ne=1.36e7,record.ftny=FALSE){
  ##browser()
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #column names
  ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  while(t_now < t){ #when current time is less than t
    rates<- rate_move(protein,protein_op,s,DisMat,MuMat,m,C,Phi,q,Ne) #moving rates of a protein to its neighbors
    lambda <- sum(rates) #total rates of moving (reaction)
    ## consider modifying this when running parallel, otherwise it'll return the same results
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
  if(record.ftny){
    ftny_vec <- NULL
    ##fitness_vec <- NULL
    ftny_vec <- sapply(1:dim(path)[1],function(i){Ftny_protein(path[i,seq(1:l)],protein_op,s,DisMat)})
    ##fitness_vec <- exp(-q*Phi*C/ftny_vec)
    path <- cbind(path,ftny_vec)
    ##path <- cbind(path,fitness_vec)
    colnames(path)[l+3] <- "Functionality"
    ##colnames(path)[l+4] <- "Fitness"
  }
  ##shift the third column up one step so that the waiting time is the time
  ## spent in the state in the first column
  path[,3] <- c(path[-1,3],NA)  
  path
}
#simulation(c(1,2),c(3,4),1000,10,0.1,GM,mumat,indep=T)


##Simulation of protein sequences of length "l" on a phylogeny "tree", 
##given the ancestral sequence "rootseq", optimal amino acid sequence "protein_op",
##selection coefficient "s", here we consider s THE SAME ACROSS ALL THE SITES IN ONE GENE 
simTree <- function(tree,l,protein_op=rep(1,l),m,s,GTRvec,alpha=al,
                    beta=be, gamma=ga,scale=T,C=2, 
                    Phi=0.5,q=4e-7, Ne=1.36e7,
                    bf=NULL,rootseq=NULL,ancestral=FALSE){
  ## check if the input tree is rooted binary tree. if not, throw an error.
  ##what if we want simulation on a branch only, i.e., there is only one tip?

  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  if(is.null(bf)) bf = rep(1/m,m) #base frequency, randomly chosen from all states, used to choose root sequence when it's not specified
  if(is.null(rootseq))rootseq = sample(c(1:m), l, replace=TRUE, prob=bf) #sequence at the root
  GM1=GM_cpv(GMcpv,alpha,beta,gamma) #Distance matrix from new weights
  mumat = aa_MuMat_form(GTRvec,scale) #mutation rate matrix, from the GTR matrix
  tree = ape:::reorder.phylo(tree) #oder the tree, cladwise
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree
  res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes
  parent <- as.integer(edge[, 1]) #parents of the edges
  child <- as.integer(edge[, 2]) #children of the edges
  root <- as.integer(parent[!match(parent, child, 0)][1]) #root node
  res[root,] <- rootseq
  tl = tree$edge.length #lengths of the edges
  ##browser()
  for(i in 1:length(tl)){
    from = parent[i] 
    to = child[i]
    ##simulation on this one branch
    for(j in 1:l){
      res[to,j] = as.numeric(tail(simulation(res[from,j],protein_op[j],tl[i],m,s,GM1,mumat,indep,
           formula,model,a1,a2,Phi,q,Ne,record.fitness=FALSE),1)[1])
      
    }
  }
  k = length(tree$tip)
  label = c(tree$tip, as.character((k+1):nNodes))
  rownames(res)=label 
  if(!ancestral)res = res[ tree$tip, , drop=FALSE]
  return(as.data.frame(res))
  ##if(pt=="AA") return(phyDat.AA(as.data.frame(res), return.index=TRUE))
}
#simTree(tree,5,protein_op=Protein_op[1:5],m=20,s=0.3,Nu_vec,al,be,ga,q=4e-7,Ne=1.36e7,Root[1:5])
