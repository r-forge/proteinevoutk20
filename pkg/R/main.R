scale <- function(Q,bf){
  l = length(bf)
  res = Q
  diag(res) = rep(0,l)
  res = res*bf
  res2 = res*rep(bf,each=l)
  diag(res) = -colSums(res)
  return (res/sum(res2))
}
ll_site <- function(tree,data,bf,Q){
  ##If the given tree is not rooted and binary, then throw error and exit
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  ##If the amino acid at the root is specified, set the base frequencies
  m = length(bf)
  Q = scale(Q,bf)
  tree <- ape:::reorder.phylo(tree,"p") #reorder the tree in pruningwise order
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree
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
  for(i in 1:tree$Nnode){ #for each interior node calculate the probability vector of observing 1 of 20 states
    from = parent[2*i] #parents
    to = child[(2*i-1):(2*i)] #direct descendents
    t_left = tl[2*i-1] #left branch length
    t_right = tl[2*i] #right branch length
    v.left <- expm(Q*t_left) #probabilities of transition from one state to another after time t
    v.right <- expm(Q*t_right)
    probvec[from,] <- as.vector((v.left%*%probvec[to[1],])*(v.right%*%probvec[to[2],])) #pruning, vector form
    check.sum <- sum(probvec[from,])
    if(check.sum==0) #probability is very very low
      warning("numerical overflow",immediate.=TRUE)
  }
  return(as.numeric(probvec[root,] %*% bf))
  #return(list(ll=max(probvec[root,]),root=which.max(probvec[root,]))) #with the corresponding root returned 
  #return(max(probvec[root,])) #just the value
}