library("qgraph")
library("network")
##set all positive entries in mumat to 1, and others to 0, save it in mumat01
net <- network(MUMAT,loops=F)
mumat01 <- net[,]

#view the color transition in stripes
Col <- rev(heat.colors(20))
image(1:20,1,as.matrix(1:20),col=Col)
##############################################################
##for a network object, find the step 1 neighbors and step 2 neighbors of node k
##return as a list with 2 entries: ng1, ng2
get.neighbor.step <- function(net,k){
  ng1 <- get.neighborhood(net,k) # 1-step neighbors of node k
  ng2 <- NULL # 2-step neighbors
  for(i in ng1){
    ngi <- get.neighborhood(net,i) # 1-step neighbors of neighbors of k
    ng2_add <- setdiff(ngi,c(ng1,k,ng2)) # 2-step neighbors of k only
    ng2 <- c(ng2,ng2_add) # put 2-step neighbors together
  }
  return(list(step1=ng1,step2=ng2))
}
##############################################################
##Given an integer, find the coordinates of that many points
spacing <- function(k){
  if(k%%2!=0){
    return((-(k-1)/2):((k-1)/2))
  }
  if(k%%2==0){
    return(setdiff((-k/2):(k/2),0))
  }
}
##Create a layout for graph according to the neighborhood status of nodes in a network
create.layout <- function(op,ng1,ng2){
  layout <- matrix(0,20,2)
  layout[op,] <- c(0,1)
  layout[ng1,] <- cbind(spacing(length(ng1)),0.5)
  layout[ng2,] <- cbind(spacing(length(ng2)),-0.5)
  return(layout)
}
##############################################################
##build a matrix with entries as difference between functionalities for corresponding nodes
##to be used for graphing
##if one state/node i is better than another j in functionality, then there is an arrow going from j to i
##i.e. if fmatrix[j,i] > 0, keep the value, otherwise, fmatrix[j,i] = 0
fmatrix <- function(ftny){
  res <- matrix(0,20,20)
  for(i in 1:20){
    for(j in 1:20){
      res[i,j] <- max(ftny[j] - ftny[i],0)
    }
  }
  return(res)
}
##############################################################
##functionality qgraph
## parameters: op - optimal amino acid; s - selection coefficient
fgraph <- function(op,s,beta=be,gamma=ga,graph=TRUE){
  ftny <- vector("numeric",length=20)
  GM <- GM_cpv(GM_CPV,al,beta,gamma)
  #functionality vector
  for(i in 1:20){
    ftny[i] <- Ftny_protein(i,op,s,GM)
  }
  ra <- rank(ftny) #rank of the functionality, used later for coloring of nodes
  f.matrix <- fmatrix(ftny) #matrix with differences between functionalities
  ng <- get.neighbor.step(net,op) # find the 1- and 2- step neighbors
  ng1 <- ng$step1
  ng2 <- ng$step2
  lay <- create.layout(op,ng1,ng2) # layout of the figure
  ##find out the out degrees (number of outflows) from each amino acid
  #the matrix that's plotted, there is edge if it's possible to mutate 
  #and only goes to higher functionality
  fm.matrix <- f.matrix*mumat01
  #number of arrows going out of the states
  # if 0 then it'll be a local optimum
  out.num <- vector("numeric",20) 
  for(i in 1:20){
  out.num[i] <- length(which(fm.matrix[i,]!=0))
  }
  num.localop <- sum(out.num==0)
  if(graph){
    ##If the outflow for a node is 0, use square for the node, otherwise circle
    node.shape <- vector(mode="character",20)
    node.shape <- ifelse(out.num==0,"square","circle")
    qgraph(f.matrix*mumat01,color=Col[ra],curve=0.3,cut=0.001,
           layout=lay,asize=0.15,shape=node.shape)
           
           ##filetype="pdf",filename=paste("A_op",op,sep=""))
  }
  return(list(out=out.num,num.localop=num.localop))
  
}
##############################################################
## given beta, gamma, and s, find the number of local optima for all optimal amino acids
## return a vector of length 20
num.localop <- function(beta,gamma,s=0.1){
  res <- lapply(1:20,fgraph,s,beta=beta,gamma=gamma,graph=F)
  sum(sapply(1:20,function(x) return(res[[x]]$num.localop>1)))
}
##vectorized version of num.localop
num.localop.vec <- Vectorize(num.localop,c("beta","gamma"))


check.local <- function(beta,gamma){
  lapply(1:20,fgraph,s=0.1,beta=beta,gamma=gamma,graph=TRUE)
}

##run function on all combinations of beta.vec and gamma.vec
# beta.vec <- seq(0,0.1,0.001)
# gamma.vec <- seq(0,0.1,0.001)
# localnum_array <- outer(beta.vec,gamma.vec,num.localop.vec)
###qgraph(fma*mumat01,color=Col[ra],curve=0.3,cut=0.0001,layout=lay,asize=0.1)
