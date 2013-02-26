# ########################################################################
# ##   GO information plots, see if attributes in the same group cluster 
# ##   together
# ########################################################################
# gos <- lapply(1:length(GGList), function(x) GGList[[x]]$order)
# sp <- scatterplot3d(sw_max[1,gos[[5]]],sw_max[2,gos[[5]]],log(sw_max[3,gos[[5]]]),
#                     xlab="g", ylab="beta",zlab="log(gamma)")
# sp$points3d(sw_max[1,gos[[2]]],sw_max[2,gos[[2]]],log(sw_max[3,gos[[2]]]),col="blue")
# sp$points3d(sw_max[1,gos[[3]]],sw_max[2,gos[[3]]],log(sw_max[3,gos[[3]]]),col="red")
# sp$points3d(sw_max[1,gos[[4]]],sw_max[2,gos[[4]]],log(sw_max[3,gos[[4]]]),col="green")
# 
# sp1 <- scatterplot3d(sw_maj[1,gos[[1]]],sw_maj[2,gos[[1]]],log(sw_maj[3,gos[[1]]]),
#                      xlab="g", ylab="beta",zlab="log(gamma)")
# sp1$points3d(sw_maj[1,gos[[2]]],sw_maj[2,gos[[2]]],log(sw_maj[3,gos[[2]]]),col="blue")
# sp1$points3d(sw_maj[1,gos[[3]]],sw_maj[2,gos[[3]]],log(sw_maj[3,gos[[3]]]),col="red")
# sp1$points3d(sw_maj[1,gos[[4]]],sw_maj[2,gos[[4]]],log(sw_maj[3,gos[[4]]]),col="green")
# 
# 
# var_s_max <- sapply(1:length(gos), function(x) var(s_max[gos[[x]]]))
# var_s_maj <- sapply(1:length(gos), function(x) var(s_maj[gos[[x]]]))
# 
# #for each gene, find the amino acids that have nonzero weights of being optimal
# # using opw_opw1 (opw_opw has one gene 25 that doesn't have a result)
# opAA <- vector(mode="list",length=106)
# for(i in 1:106){
#   index <- which(opw_opw1[,i]!=0)
#   value <- opw_opw1[index,i]
#   opAA[[i]]$AA <- index[order(value,decreasing=TRUE)]
#   opAA[[i]]$weights <- value[order(value,decreasing=TRUE)]
# }
# 
# s_list <- vector(mode="list",length=length(gos))
# for(i in 1:length(gos)){
#   s_list[[i]] <- s_max[gos[[i]]]
# }
# beta_list <- vector(mode="list",length=length(gos))
# for(i in 1:length(gos)){
#   beta_list[[i]] <- GM_max[2,gos[[i]]]
# }
# gamma_list <- vector(mode="list",length=length(gos))
# for(i in 1:length(gos)){
#   gamma_list[[i]] <- GM_max[3,gos[[i]]]
# }
# 
# s <- 1.486043
# beta <- 0.1327872 
# gamma <- 0.0009136954
# Q <- c(4.028232, 18.44403, 10.54094, 19.01206, 4.451557, 1)
# brlen <- c(0.2590216, 0.126402, 0.03749147, 0.5138506, 0.1755939, 0.4282816, 0.09928871, 0.2314689, 6.116373, 2.87642, 8.022745, 0.8643566)
# 
# sitelikmat <- matrix(0,nrow=107,ncol=20)
# for(i in 1:20){
#   rootbf <- rep(0,20)
#   rootbf[i] <- 1
#   resi <- mllm(data,treeir,s,beta,gamma,Q,opaa=res$ll$opaa,rootbf=rootbf)
#   sitelikes <- resi$ll$sitelik
#   sitelikmat[,i] <- sitelikes
# }
# siteprob <- exp(sitelikmat)
# siteprob = apply(siteprob,1,function(x) x/sum(x))
# indexp = attr(data,"index")
# siteprob = siteprob[,indexp]
# 
# startsq <- apply(siteprob,2,sample,x=20,size=1,replace=TRUE)

######################################################################
# ## Read the WAG matrix (lower triangular part)
# WagMat <- scan("~/proteinevoutk20/pkg/Data/wag.txt")
# #Qwag <- mat_form_lowtriQ(Q=WagMat,bf=bfaa)
# ll_site_lowQ <- function(tree,data,Q,bf=rep(1/20,20),g=1){
#   ##If the given tree is not rooted and binary, then throw error and exit
#   if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
#   m = 20
#   ##if the base frequencies are not specified, do a uniform distribution
#   if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
#   Qmat <- mat_form_lowtriQ(Q,bf)
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
#   #return(list(ll=max(probvec[root,]),root=which.max(probvec[root,]))) #with the corresponding root returned 
#   #return(max(probvec[root,])) #just the value
# }
# 
# simulationQ <- function(protein,t,Q,bf=rep(1/20,20)){
#   Qmat <- mat_form_lowtriQ(Q,bf)
#   l <- length(protein) #number of sites
#   t_now <- 0 #time until the current step of simulation
#   path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
#   colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #column names
#   ##last two columns in the paths, recording Time up to this point and the waiting time at the current state
#   colnames(path)[l+1] <- "Time Now"
#   colnames(path)[l+2] <- "Waiting Time"
#   m = 20
#   while(t_now < t){ #when current time is less than t
#     rates <- Qmat[protein,]
#     rates <- c(t(rates)) #moving rates of a protein to its neighbors
#     rates[rates < 0] <- 0
#     lambda <- sum(rates) #total rates of moving (reaction)
#     t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
#     t_now <- t_now + t_wait # current time after the first transition
#     ##determine where the first transition is
#     if(t_now < t){
#       index <- sample(1:length(rates),1,replace=T,rates) #index of the protein it moves to
#       pos <- (index-1) %/% m + 1 #index of the site that changes
#       rmd <- (index-1) %% m + 1 # the amino acid that site changes to
#       protein[pos] <- seq(1:m)[rmd] #new protein
#       path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
#     }
#     else{
#       path <- rbind(path,c(protein,t,NA)) #record this moving step
#     }
#   }
#   ##shift the third column up one step so that the waiting time is the time
#   ## spent in the state in the first column
#   path[,l+2] <- c(path[-1,l+2],NA)  
#   return(path)
# }
# 
# ## simulation based on WAG model, and the plots of functionality along the path using the simulations
# sim.Wag <- function(t=brlen,protein,bf=bfaa){
#   sim <- simulationQ(protein=protein,t=t,Q=WagMat,bf=bf)
#   l <- dim(sim)[2] - 2
#   fty <- apply(sim[,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
#   dis <- apply(sim[,1:l],MARGIN=1,FUN=pchem_d,protein2=opaa,DisMat=dismat)#distance from optimal amino acids
#   dis <- apply(dis,2,mean)#average distance for all sites
#   ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
#   disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
#   #plot(ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
#   return(list(sim=sim,fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
# }
# 
# sim.New <- function(s=1,t=10,root=root,opaa=opaa,beta,gamma,bfaa=bfaa){
#   dismat = GM_cpv(GM_CPV,al,beta,gamma)
#   mumat = aa_MuMat_form(res_op$Q)
#   sim <- simulation(root,opaa,t=t,s=s,DisMat=dismat,MuMat=mumat,bfaa=bfaa) #simulation
#   l <- dim(sim)[2] - 2
#   fty <- apply(sim[,1:l],MARGIN=1,FUN=Ftny_protein,protein_op=opaa,s=s,DisMat=dismat)#functionality
#   dis <- apply(sim[,1:l],MARGIN=1,FUN=pchem_d,protein2=opaa,DisMat=dismat)#distance from optimal amino acids
#   dis <- apply(dis,2,mean)#average distance for all sites
#   ftyfun <- stepfun(sim[-1,l+1],fty,f=0,right=FALSE)#make step functions
#   disfun <- stepfun(sim[-1,l+1],dis,f=0,right=FALSE)
#   #if(func) #plot functionality vs time
#     #plot(ftyfun,xlab="time",ylab="functionality",main=paste("functionality, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
#   #if(dist) #plot distance vs time
#     #plot(disfun,xlab="time",ylab="distance",main=paste("distance, s=",round(s,3),sep=""),pch=20,xlim=c(0,t),xaxs="i",add=add)
#   #return(as.numeric(tail(sim,1)[1:l])) #the sequence at the end of simulation
#   return(list(sim=sim,fty=fty,dis=dis,ftyfun=ftyfun,disfun=disfun)) #store the simulation result for later use
# }
######################################################################
maxdir <- "~/BackupProEvo/Newton/rokas_max/"
res_max <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maxdir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_max[[genect]] <- res_op
}

source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/hessian.R")

std <- matrix(0,nrow=3,ncol=l)
for(genect in 1:l){
  res_op <- res_max[[genect]]
  hes <- find_hessian(c(res_op$s,res_op$GMweights[2:3]),res_op$data,res_op$tree,Q=res_op$Q)
  std[,genect] <- sqrt(diag(solve(hes)))
}
save(std,file="~/proteinevoutk20/pkg/Result/std.RData")