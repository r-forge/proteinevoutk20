# rm(list=ls())
# gene = 20
# RDatafile <- paste("~/BackupProEvo/Newton/rokas7/simAna/Eqm/gene",gene,".RData",sep="")
# load(RDatafile)
# source("~/proteinevoutk20/pkg/R/prune.R")
# source("~/proteinevoutk20/pkg/R/plot_simtree.R")
# ftny.vec <- apply(datanum,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
# ## find the starting and ending position of all edges in the tree
# #plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=TRUE,plotdis=FALSE)
# plot_trace(sim=sim,plottip=TRUE,plotftny=TRUE,plotdis=FALSE)
# #plot_trace(sim=sim$emp[[1]],plottip=FALSE,plotftny=FALSE,plotdis=TRUE)
# plot_trace(sim=sim,plottip=TRUE,plotftny=FALSE,plotdis=TRUE)
#################################################
# s = 5
# fixmatall <- fixmatAll(s,DisMat=dismat)
# Qall = QAllaa1(fixmatall,mumat)
# bf <- sapply(Qall,FUN=eqmQ)
# for(i in 1:20)
#   plot(sort(bf[,i],decreasing=TRUE),type="s",ylim=c(0,1))
#################################################
load("~/proteinevoutk20/pkg/Data/Rokas/rokas.RData")
data = rokas[,seq(AArange[1,1],AArange[1,2])]
dataphy <- phyDat(data,type="AA")
tree <- read.nexus("~/proteinevoutk20/pkg/Data/GTR.tre")
s <- 1
dismat <- GM_cpv(GM_CPV,al,be,ga)
fixmatall <- fixmatAll(s,DisMat=dismat)
Q <- rep(1,6)
mumat = aa_MuMat_form(Q)
Qall <- QAllaa1(fixmatall,mumat)
ll3m.test <- function (dat1, tree, scale = 1,ancestral = NULL, ancStates = NULL,Q, g = 1) 
{
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorderPruning(tree)
  q = length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m = length(edge) + 1
  Q = Q*scale
  Q = t(Q) # in order to use C function in phangorn package
  el <- tree$edge.length
  P <- getPm(el,Q,g)
  nr <- as.integer(attr(dat1, "nr"))
  nc <- as.integer(attr(dat1, "nc"))
  node = as.integer(node - min(node))
  edge = as.integer(edge - 1)
  nTips = as.integer(length(tree$tip))
  mNodes = as.integer(max(node) + 1)
  contrast = attr(dat1, "contrast")
  nco = as.integer(dim(contrast)[1])
  res <- .Call("LogLik2", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
#   if((!is.null(ancestral))&&(!is.null(ancStates)))
#     warning("both ancestral state frequencies and states are specified, only ancestral states are used!")
## if root state frequencies and root states are both given, root states have priority and are used
  #browser()
  if(!is.null(ancStates)){
    if(length(ancStates)!=nr) #length has to be the same as the dinstince site patterns
      stop("ancestral state sequence has wrong length!")
    result = sapply(1:nr,function(x) res[[1]][x,ancStates[x]])
    result = matrix(log(result),ncol=1) #take logarithm 
    ancStates = ancStates
  }
  else{
    if(is.null(ancestral)){
      result = sapply(1:nr,function(x) max(res[[1]][x,])) #find the root state that maximizes the likelihood
      result = matrix(log(result),ncol=1) #take logarithm 
      ancStates =  sapply(1:nr,function(x) which.max(res[[1]][x,])) # return the root states at all sites
    }
    else{
      ancestral = as.matrix(ancestral,nrow=20)
      if(ncol(ancestral)==1)
        ancestral = sapply(1:nr,function(x) ancestral)
      ancestral <- apply(ancestral,2,function(x) x/sum(x))
      result = log(diag(res[[1]]%*%ancestral)) #updated phangorn, phangorn 1.7
    }
  }
  if(is.null(ancestral)){
    diag.mat = diag(20)
    ancestral = diag.mat[,ancStates]
  }
  return(list(result=result,ancestral=ancestral,ancStates=ancStates))
  #return(result)
}
