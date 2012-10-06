## source the R file that contains required functions
#source("~/proteinevoutk20/pkg/R/pphyproevo.R")
########################################################################################################
## Given s, Distance matrix, mutation rate matrix, and other parameters, find the list of 20 rate matrices corresponding to optimal aa 1:20
## same function as QAllaa, but all the matrices are scaled in this function
rate_move_mat <- function(s, DisMat, MuMat,bfaa = rep(1/20,20),C=2, Phi=0.5, q=4e-7, Ne=1.36e07){
  fn = function(i){
    mat = mat_gen_indep(i,s,DisMat,MuMat,C,Phi,q,Ne)
    mat = scaleQ(mat,bfaa)
  }
  mat_list = lapply(1:20,fn)
  mat_list
}
rate_move <- function(protein, protein_op,s, DisMat, MuMat, bfaa = rep(1/20,20),C=2, Phi=0.5, q=4e-7, Ne=1.36e07){
    mat = rate_move_mat(s,DisMat,MuMat,bfaa,C,Phi,q,Ne)
    vec_list = sapply(1:length(protein),function(i) {mat[[protein_op[i]]][protein[i],]},simplify="array")
    vec_list[vec_list < 0] = 0
    return(c(vec_list))
}
########################################################################################################
#Simulation. Given the starting protein and the optimal protein
#t should be equal to expected number of substitutions
simulation <- function(protein,protein_op,t,s,DisMat,MuMat,bfaa=rep(1/20,20),C=2, Phi=0.5,q=4e-7, Ne=1.36e7,record.ftny=FALSE){
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #column names
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
  if(record.ftny){
    ftny_vec <- NULL
    ftny_vec <- sapply(1:dim(path)[1],function(i){Ftny_protein(path[i,seq(1:l)],protein_op,s,DisMat)})
    path <- cbind(path,ftny_vec)
    colnames(path)[l+3] <- "Functionality"
  }
  ##shift the third column up one step so that the waiting time is the time
  ## spent in the state in the first column
  path[,l+2] <- c(path[-1,l+2],NA)  
  path
}
#simulation(c(1,2),c(3,4),1000,10,0.1,GM,mumat,indep=T)
#mllm(data,tree,s,beta,gamma,Q=NULL,opw=NULL,bfnu=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7)

##Simulation of protein sequences of length "l" on a phylogeny "tree", 
##given the ancestral sequence "rootseq", optimal amino acid sequence "protein_op",
##selection coefficient "s", here we consider s THE SAME ACROSS ALL THE SITES IN ONE GENE 
simTree <- function(tree,l,protein_op,s,GTRvec,alpha=al,beta=be, gamma=ga,bfaa=NULL,
                    C=2,Phi=0.5,q=4e-7, Ne=1.36e7,rootseq=NULL,ancestral=FALSE,simple=FALSE){
  call = match.call()
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  if(is.null(bfaa)) bfaa = rep(1/20,20) #base frequency, randomly chosen from all states, used to choose root sequence when it's not specified
  m=20 # number of amino acids
  if(is.null(rootseq)) rootseq = sample(c(1:m), l, replace=TRUE, prob=bfaa) #sequence at the root
  GM1=GM_cpv(GMcpv,alpha,beta,gamma) #Distance matrix from new weights
  mumat = aa_MuMat_form(GTRvec) #mutation rate matrix, from the GTR matrix
  mumat = sym.to.Q(mumat,bfaa)
  if (is.null(attr(tree, "order")) || attr(tree, "order") !="cladewise") 
    tree <- ape:::reorder.phylo(tree)
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree
  res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes
  parent <- as.integer(edge[, 1]) #parents of the edges
  child <- as.integer(edge[, 2]) #children of the edges
  root <- as.integer(parent[!match(parent, child, 0)][1]) #root node
  res[root,] <- rootseq
  tl = tree$edge.length #lengths of the edges
  for(i in 1:length(tl)){
    from = parent[i] 
    to = child[i]
    ##simulation on this one branch
    if(simple){
      for(j in 1:l){
        res[to,j] = as.numeric(tail(simulation(res[from,j],protein_op[j],tl[i],s,GM1,mumat,bfaa,
             C,Phi,q,Ne,record.ftny=FALSE),1)[1])
      }
    }
    else{
      res[to,] = as.numeric(tail(simulation(res[from,],protein_op,tl[i],s,GM1,mumat,bfaa,
                        C,Phi,q,Ne,record.ftny=FALSE),1)[1:l])
    }
  }
  res = AA[res]
  res = matrix(res,ncol=l)
  k = length(tree$tip)
  label = c(tree$tip, as.character((k+1):nNodes))
  rownames(res)=label 
  if(!ancestral)res = res[ tree$tip, , drop=FALSE]
  #res=as.data.frame(res)
  return(list(data=res,tree=tree,optimal=protein_op,rootseq=rootseq,s=s,Q=GTRvec,GMweights=c(alpha,beta,gamma),bfaa=bfaa,call=call))
  ##if(pt=="AA") return(phyDat.AA(as.data.frame(res), return.index=TRUE))
}
#simTree(tree,5,protein_op=Protein_op[1:5],m=20,s=0.3,Nu_vec,al,be,ga,q=4e-7,Ne=1.36e7,Root[1:5])
