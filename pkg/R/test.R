nu_comb <- combinations(4,2,v=Nu) #combinations of nucleotides
nu_rates <- apply(nu_comb,MARGIN=1,c2s) #convert numeric pairs into characters "ac","ag", etc
# 
###########################################################################
## Matrix of dimension 61 by 61, each entry is the index of mutation rate 
## between the 2 codons in the 6-element vector NU_rates
diff_pos <- matrix(0,61,61)
diag(diff_pos) <- 0
for(i in 1:60){
  for(j in (i+1):61){
    cmp <- cdcmp(s2c(CDS[i]),s2c(CDS[j]))
    if(cmp$num==1){
      nu_from <- s2c(CDS[i])[cmp$pos]
      nu_to <- s2c(CDS[j])[cmp$pos]
      nu_pair <- c2s(sort(c(nu_from,nu_to)))
      diff_pos[i,j] <- as.numeric(factor(nu_pair,levels=nu_rates,labels=1:6))
    }
    else
      diff_pos[i,j] <- 0
  }
}
diff_pos <- diff_pos + t(diff_pos)
dimnames(diff_pos)[[1]] <- CDS #rownames and colnames are the amino acid names
dimnames(diff_pos)[[2]] <- CDS
###########################################################################
## indAApair: a list of length 190, each corresponds to a pair of amino acids
## including the number of possible ways of chaning to each other, and the rates for all the ways
AApairs <- combinations(20,2)
AApairs <- sapply(1:190,function(x) AApairs[x,],simplify=FALSE)
indAApair <- vector("list",length(AApairs))
for(i in 1:length(AApairs)){
  diffmat <- diff_pos[cdlist[[AApairs[[i]][1]]],cdlist[[AApairs[[i]][2]]]]
  if(sum(diffmat)!=0){
    indAApair[[i]]$rates <- as.numeric(diffmat[diffmat!=0]) #index of rates in NU_rates
    indAApair[[i]]$tot <- length(diffmat) #total number of codon comb for the pair of AA's
  }
  else{
    indAApair[[i]]$rates <- 0
    indAApair[[i]]$tot <- 1
  }
  indAApair[[i]]$AA <- AA[AApairs[[i]]]
}

###########################################################################
## look at the nonzero entries in the AArate matrix, and try to figure out the 
## influence of NU_rates on AA_rates
nonzeroPair <- indAApair[!sapply(indAApair,FUN=is.null)]
sapply(1:length(nonzeroPair),function(x) nonzeroPair[[x]]$rates)
sapply(1:length(nonzeroPair),function(x) nonzeroPair[[x]]$tot)
NU_to_AA <- matrix(0,nrow=6,ncol=length(nonzeroPair))
for(i in 1:length(nonzeroPair)){
  rates <- nonzeroPair[[i]]$rates
  NU_to_AA[,i] <- sapply(1:6,function(x) sum(rates==x))/nonzeroPair[[i]]$tot
}


# get.root <- function(mat){
#   l <- dim(mat)[2]
#   sapply(1:l,function(x) which(mat[,x]==1))
# }


# simulationQ1 <- function(protein,t,Q=NULL,bf=NULL){
#   if(is.null(bf)) bf <- rep(1/20,20)
#   if(is.null(Q)) Q <- rep(1,190)
#   if(is.matrix(Q)) Qmat <- Q
#   else if(is.vector(Q))  
#     Qmat <- mat_form_lowtriQ(Q,bf,byrow=TRUE) #transition rate matrix
#   l <- length(protein) #number of sites
#   P <- expm(Qmat*t)
#   res <- protein
#   for(i in 1:20){
#     index <- which(protein==i)
#     res[index] <- sample(1:20,length(index),replace=TRUE,prob=P[i,])
#   }
#   return(res)
# }
# 
# find.Q <- function(seq1,seq2,Q=rep(1,6),t=1,method="COBYLA",maxeval="1000",print_level=0){
#   fn <- function(Q){
#     Q <- exp(Q)
#     mumat <- aa_MuMat_form(Q)
#     P <- expm(mumat*t)
#     prob <- -log(sapply(1:length(seq1), function(i) P[seq1[i],seq2[i]]))
#     cat(Q,sum(prob),"\n")
#     return(sum(prob))
#   }
#   lower=rep(-Inf,6)
#   upper=rep(Inf,6)
#   #options for optimizer
#   opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
#                "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
#   res = nloptr(x0=log(Q),eval_f=fn, lb=lower,ub=upper,opts=opts)
#   #optimx(log(Q),fn,hessian=FALSE,method=method,itnmax=itnmax)
#   return(res)
# }
# 
# find.Q1 <- function(seq1,seq2,q=1,t=1,interval=c(-10,10)){
#   fn <- function(q){
#     q = exp(q)
#     Q <- rep(q,6)
#     mumat <- aa_MuMat_form(Q)
#     P <- expm(mumat*t)
#     prob <- -log(sapply(1:length(seq1), function(i) P[seq1[i],seq2[i]]))
#     cat(q,sum(prob),"\n")
#     return(sum(prob))
#   }
#   optimize(f=fn,interval=interval)
# }
# 
# optim.ultrametric <- function(tree,data,method="SBPLX",maxeval="100",print_level=1,...){
#   ## free parameters: total tree depth, all the internal branch lengths
#   ## calculate the exterior branch lengths given the free parameters
#   nTips <- length(tree$tip.label)
#   nIntBr <- nTips - 2 #number of interior branches
#   edge <- tree$edge
#   int.ind <- which(!(edge[,2] %in% c(1:nTips))) #index of interior branches in tree$edge
#   ext.ind <- which(edge[,2] %in% c(1:nTips))
#   fn <- function(x){ ## x <- c(intBrlen, depth)
#     x <- exp(x)
#     cat(x,"\n")
#     depth <- x[nIntBr + 1] # tree depth is the last element in the vector
#     tree$edge.length[int.ind] <- x[1:nIntBr] # assign internal branches new lengths
#     tree <- ult.tree(tree,depth=depth)$tree #calculate exterior branch lengths
#     cat("tree edge lengths","\n")
#     cat(tree$edge.length,"\n")
#     if(any(tree$edge.length<0)) return(10^7)
#     else{
#       res <- -mllm1(data,tree,...)$ll$loglik
#       cat("loglik = ", res, "\n")
#       return(res)
#     }
#   }
#   ## constraint function: add "eval_g_ineq=g0" in the nloptr function
# #   g0 <- function(x){
# #     x <- exp(x)
# #     depth <- x[nIntBr + 1]
# #     tree$edge.length[int.ind] <- x[1:nIntBr] # assign internal branches new lengths
# #     tree <- ult.tree(tree,depth=depth) #calculate exterior branch lengths
# #     return(-tree$edge.length[ext.ind])
# #   }
#   lower=rep(-Inf,nIntBr+1)
#   upper=rep(Inf,nIntBr+1)
#   #options for optimizer
#   opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
#                "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
#   
#   res = nloptr(x0=log(c(tree$edge.length[int.ind],max(node.depth.edgelength(tree)))),
#                eval_f=fn, lb=lower,ub=upper,opts=opts)
#   x <- res$solution
#   x <- exp(x)
#   depth <- x[nIntBr + 1] # tree depth is the last element in the vector
#   tree$edge.length[int.ind] <- x[1:nIntBr] # assign internal branches new lengths
#   tree <- ult.tree(tree,depth=depth)$tree #calculate exterior branch lengths
#   return(list(res=res,tree=tree))
# }

# internal.br <- function(tree){
#   nTips <- length(tree$tip.label)
#   edge <- tree$edge
#   internal.ind <- which(!(edge[,2] %in% c(1:nTips)))
#   return(list(index = internal.ind,edge=edge[internal.ind,],edge.length=tree$edge.length[internal.ind]))
# }


