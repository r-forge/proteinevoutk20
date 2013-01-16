library("qgraph")
library("network")
source("~/proteinevoutk20/pkg/R/main.R")
Col <- rev(heat.colors(20))
##############################################################
# ##A test plot using "network" package and/or "qgraph" package## --not used later, just a test
# A <- mat_gen_indep(1,0.1,GM,MUMAT) #transition rate matrix
# #sides <- as.numeric(table(aa)) # number of sides equals number of codons for a particular aa
# ftny <- sapply(1:20,Ftny_protein,protein_op=1,s=0.1,GM)
# plot(net,edge.lwd=A*10,displaylabels=T,vertex.cex=8-log(ftny,0.1)*100)
# #qgraph(A1,color=sides+1,vsize=(7.75-log(ftny,0.1)*100)*0.8,esize=8,zsize=0.1,filetype="pdf")
# qgraph(A,vsize=(7.75-log(ftny,0.1)*100)*0.8,esize=8,zsize=0.1)
# # 
# ##############################################################
##qgraph the substitution rate matrix, given the optimal aa, s, and expoent in Ne and the desired layout
##color represents the funtionality
##return the qgraph object, with information on layout and etc.
aa_plot <- function(op=1,s=1,beta=be,gamma=ga,bfaa=rep(1/20,20),Q=NU_VEC,layout="spring"){
  GM <- GM_cpv(GM_CPV,alpha=al,beta=beta,gamma=gamma)
  mumat = aa_MuMat_form(Q)
  ftny <- sapply(1:20,Ftny_protein,protein_op=op,s=s,DisMat=GM)
  ra <- rank(ftny)
  #plot(net,edge.lwd=A*10,displaylabels=T,vertex.col=sides,vertex.cex=8-log(ftny,0.1)*100)
  A <- mat_gen_indep(op,s,GM,mumat)
  A <- scaleQ(A,bf=bfaa)
  gpara <- qgraph(A,color=Col[ra],esize=10,zsize=0.1,layout=layout,edge.color="black",cut=0.06,
                  details=TRUE,labels=AA)
                  #filetype="pdf",filename=paste("A_op",op,"_Ne",exponent,"_s",s,sep=""))
  return(gpara)
}
aa_plot(1,10,beta,gamma,bfaa,res_op$Q,layout=qplot$layout)

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
