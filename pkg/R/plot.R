library("qgraph")
library("network")
##############################################################
# ##A test plot using "network" package and/or "qgraph" package## --not used later, just a test
# A <- mat_gen_indep(1,0.1,GM,mumat,q=4e-7,Ne=1.36e7) #transition rate matrix
# sides <- as.numeric(table(aa)) # number of sides equals number of codons for a particular aa
# ftny <- vector("numeric",length=20) #functionality vector
# op <- 1 #optimal aa
# s <- 0.1 #selection coefficient
# for(i in 1:20){
#   ftny[i] <- Ftny_protein(i,op,s,GM)
# }
# plot(net,edge.lwd=A*10,displaylabels=T,vertex.col=sides,vertex.cex=8-log(ftny,0.1)*100)
# #qgraph(A1,color=sides+1,vsize=(7.75-log(ftny,0.1)*100)*0.8,esize=8,zsize=0.1,filetype="pdf")
# qgraph(A,color=sides+1,vsize=(7.75-log(ftny,0.1)*100)*0.8,esize=8,zsize=0.1)
# 
# ##############################################################
# 
# 
# ##qgraph the substitution rate matrix, given the optimal aa, s, and expoent in Ne and the desired layout
# ##color represents the funtionality (green->purple, worst -> best)
# ##return the qgraph object, with information on layout and etc.
# aa_plot <- function(op,s,exponent,layout="spring"){
#   ftny <- vector("numeric",length=20)
#   #functionality vector
#   for(i in 1:20){
#     ftny[i] <- Ftny_protein(i,op,s,GM)
#   }
#   ra <- rank(ftny)
#   #plot(net,edge.lwd=A*10,displaylabels=T,vertex.col=sides,vertex.cex=8-log(ftny,0.1)*100)
#   A <- mat_gen_indep(op,s,GM,mumat,q=4e-7,Ne=1.36*10^exponent)
#   gpara <- qgraph(A,color=Col[ra],esize=8,zsize=0.1,layout=layout)
#                   #filetype="pdf",filename=paste("A_op",op,"_Ne",exponent,"_s",s,sep=""))
#   return(gpara)
# }
# ##############################################################
# #s <- 0.1
# #exponent <- 7
# 
# ###Plots###
# ##For each amino acid as the optimal one, given selection coeff and Ne
# ##Graph the subsitution rate matrix, with the sizes of nodes represent
# ##the proportion of time spending in the corresponding states from simulations
# ##Simulation is done with starting point from any of 20 amino acids
# for(op in 3:20){
#   ftny <- vector("numeric",length=20)
#   #functionality vector
#   for(i in 1:20){
#     ftny[i] <- Ftny_protein(i,op,s,GM)
#   }
#   ra <- rank(ftny) #rank of the functionality, used later for coloring of nodes
#   A <- mat_gen_indep(op,20,s,GM,mumat,q=4e-7,Ne=1.36*10^exponent) #substitution rate matrix
#   
#   ##get the layout
#   ng <- get.neighbor.step(net,op)
#   ng1 <- ng$step1
#   ng2 <- ng$step2
#   lay <- create.layout(op,ng1,ng2)
#   ##Plot the network, color of node represents the functionality, edge colors correspond to the rates
#   ##gpara <- qgraph(A,color=Col[ra],esize=8,zsize=0.1,
#                   #filetype="pdf",filename=paste("A_op",op,"_Ne",exponent,"_s",s,sep=""))
#   ##Do simulations with different starting points, and then vary the size of nodes
#   ##according to the freqencies they appear in the simulated chain
#   for(st in 1:20){
#     sim <- simulation(st,op,500,20,0.1,GM,mumat,q=4e-7,Ne=1.36*10^exponent)
#     ##frequencies of different states
#     ##freq <- getfreq(sim)
#     tfreq <- getime(sim) #This function find the amounts of time spent in all states from simulation result
#     tfreq <- tfreq/max(tfreq)
#     #ra1 <- rank(tfreq)
#     qgraph(A,color=Col[ra],esize=8,zsize=0.1,vsize=tfreq*7,layout=lay,cut=0.01,asize=0.1,
#            filetype="pdf",filename=paste("A_op",op,"_Ne",exponent,"_s",s,"_",st,sep=""))
#   }
# }
# ##############################################################
# #Do the transition plots for different optimal aminio acids, different selection forces and population sizes
# #The layout for each optimal aa is fixed at the same so that it's easy to compare them
# #layout here is "spring" layout
# #OP, S and EX are the vector where op, s and ex are drawn from, respectively
# #ex. OP <- c(1:20), s <- c(0.01,0.1,1,10), EX <- c(6:10)
# PlotSpring <- function(OP,S,EX){
#   for(op in OP){
#     a <- aa_plot(op,0.1,7)
#     for(s in S){
#       for(ex in EX){
#         aa_plot(op,s,ex,layout=a$layout)
#       }
#     }
#   }
# }
# ##############################################################
# ### find the time spent in all differnt amino acids
# getime <- function(data){
#   l <- dim(data)[1]
#   data <- data[c(-l,-l+1),]
#   timesp <- vector("numeric",20) #the vector to store the result
#   for(i in 1:20){
#     if(i %in% data[,1]){ #if ith state appears in the simulation
#       data_i <- data[data[,1]==i,] #part of the data with i as the state
#       if(length(data_i)==3){ #if this state only appears once
#         timesp[i] <- data_i[3]
#       }
#       else{ #more than once
#         timesp[i] <- sum(data_i[,3])
#       }
#     }
#     else{ #the state is not hit in the simulation
#       timesp[i] <- 0
#     }
#   }
#   return(timesp)
# }
# 
# getfreq <- function(data) 
# {
#     freq <- as.numeric(table(data[, 1]))
#     #freq <- freq/max(freq)
#     return(freq)
# }
# ##############################################################

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
