#MLE for s, given beta and gamma and all other parameter values
optim.s <- function(data, tree, s,method="BOBYQA",maxeval="100",print_level=0, ...){ # s is the initial
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,...)
  #these don't change with the change of s
  dismat = res.initial$dismat
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  fn = function(ab,data,tree){
    result = -mllm1(data=data,tree=tree,s=ab,dismat=dismat,mumat=mumat,bfaa=bfaa, ...)$ll$loglik
    return(result)
  }
  lower <- 0 #lower bound
  upper <- Inf #upper bound
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,
               "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=s,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}
# find mle of s for a range of genes in rokas data, given values of beta and gamma
## only apply to rokas data, for other data sets, need to tweak the codes
optim.s.range<- function(beta,gamma,generange,tree,multicore=FALSE,method="BOBYQA",maxeval="100",print_level=0, ...){
  mle.s.one <- function(k){ ## find mle of s for one gene
    ## for now all search start with initial value "1" for s, could let user control the starting value             
    mle <- optim.s(ROKAS_DATA[[k]],tree,s=1,method=method,maxeval=maxeval,print_level=print_level,beta=beta,gamma=gamma,...)
    return(mle)
  }
  if(multicore)
    res <- mclapply(generange,mle.s.one)
  else
    res <- lapply(generange,mle.s.one)
  
  return(res)
}
######################################################################################################
# find mle for beta and gamma, that maximize the likelihood for all genes,Q is given
optim.w <- function(beta,gamma,generange,tree,multicore=FALSE,method="BOBYQA",maxeval="100",print_level=0, ...){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  ab <- c(beta,gamma)
  fn <- function(ab,generange,tree){
    cat("beta, gamma = ", ab, "\n")
    #call the previous function to optimize s for all genes
    mle <- optim.s.range(ab[1],ab[2],generange,tree=tree,multicore=multicore,method=method,maxeval=maxeval,print_level=print_level,...) 
    mle.val <- sapply(1:length(generange),function(ind) mle[[ind]]$objective) # best -loglikelihood values
    return(sum(mle.val)) #summation of all values
  }
  
  lower <- c(0,0) #lower bound
  upper <- c(Inf,Inf) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,generange=generange,tree=tree)
  return(res)
}

######################################################################################################
#MLE for s, beta and gamma, using subplex method by default
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,Q=NU_VEC))
optim.s.weight <- function(data, tree, s,beta,gamma, method="SBPLX",maxeval="50",print_level=0,...){
  #store information from initial condition, with other parameters fixed
  if(is.null(s)) s = 1
  if(is.null(beta)) beta=be
  if(is.null(gamma)) gamma=ga
  res.initial = mllm1(data=data,tree=tree,s=s,beta=beta,gamma=gamma,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  ab <- c(s,beta,gamma) ##initial value
  #   ab[ab<=0] <- 1e-4
  #   ab <- log(ab)
  fn = function(ab,data,tree){
    #ab <- exp(ab)
    #print(ab)
    if(any(ab > 20)) return(10^5)
    else{
      result = -mllm1(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], mumat=mumat, ...)$ll$loglik
      return(result)
    }
  }
  lower <- rep(0,3)
  upper <- rep(10,3)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}
######################################################################################################
#MLE for s, beta and gamma, using subplex method by default
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,scale.vec=rep(1,20)
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,Q=NU_VEC))
optim.scale <- function(data, tree, scale.vec=rep(1,20),Qall=NULL,method="SBPLX",maxeval="50",print_level=0,...){
  if(is.null(Qall)){
    res.initial = mllm1(data=data,tree=tree,scale.vec=scale.vec,...)
    Qall <- res.initial$Qall
  }
  scale.vec <- scale.vec/scale.vec[1]
  ab <- scale.vec[-1] ##initial value
  #   ab[ab<=0] <- 1e-4
  #   ab <- log(ab)
  fn = function(ab,data,tree){
    #ab <- exp(ab)
    #print(ab)
    scale.vec <- c(1,ab)
    result = -mllm1(data=data,tree=tree,scale.vec=scale.vec,Qall=Qall,...)$ll$loglik
    return(result)
  }
  lower <- rep(0,19)
  upper <- rep(Inf,19)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}
######################################################################################################
## MLE for opw (weights for optimal amino acids)
## supply the function, the gradient function, and the restriction function (sum=1)
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call: optim.opw(data,tree,opw=rep(1,20),s=0.1,beta=beta,gamma=gamma,Q=NU_VEC...)
optim.opw <- function(data, tree,opw=NULL,print_level=0, ...){
  if(is.null(opw))
    opw = findBf2(data)
  opw = opw/sum(opw) #opw given in the function call doesn't have to sum to 1
  
  res = mllm(data=data,tree=tree,opw=opw,...) #store llmat (loglikelihood values for all opaa) 
  llmat = exp(res$ll$llmat)    #so that they don't need to be evaluated again and again
  weight = attr(data,"weight")
  #cat("opw optimization, starting loglikelihood = ", res$ll$loglik, "\n") #function value at the starting point
  
  # function to optimize on and its gradient function
  eval_f_list <- function(opw){
    return(list("objective"=llaaw1(opw,weight=weight,llmat=llmat),
                "gradient"=llaaw_grad(opw,weight=weight,llmat=llmat)))
  }
  # linear equality constraint, and its jacobian
  eval_g_list <- function(opw){
    return(list("constraints"=sum(opw)-1,"jacobian"=rep(1,20)))
  }
  lower <- rep(0,20) #lower bound
  upper <- rep(1,20) #upper bound
  local_opts <- list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1e-7) #options for local optimizer
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","maxeval"="1000000","xtol_rel"=1e-7,"ftol_rel"=.Machine$double.eps,
               "local_opts"=local_opts,"print_level"=print_level)
  res = nloptr(x0=opw,eval_f=eval_f_list,eval_g_eq=eval_g_list, lb=lower,ub=upper,opts=opts)
  # res$objective: best function value found
  # res$solution: best parameter values
  return(res)
}

# find mle for beta and gamma, that maximize the likelihood for all genes,Q is given
optim.opw.range <- function(s,beta,gamma,generange,Q,tree,multicore=FALSE,print_level=0){
  if(multicore)
    res <- mclapply(ROKAS_DATA[generange],optim.opw,tree=tree,opw=rep(1,20),
                    print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q)
  else
    res <- lapply(ROKAS_DATA[generange],optim.opw,tree=tree,opw=rep(1,20),
                  print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q)
  return(res)
}
## optimize beta, gamma and s, for each gene find the best opw for them, instead of using a universal one
optim.opw.sbg <- function(s,beta,gamma,Q,tree,generange,multicore=FALSE,print_level=0,...){
  ab <- c(s,beta,gamma)
  fn <- function(ab){
    cat("s,beta,gamma = ", ab, "\n")
    mle <- optim.opw.range(ab[1],ab[2],ab[3],generange=generange,Q=Q,tree=tree,multicore=multicore,
                           print_level=print_level,...)
    mle.val <- sapply(1:length(generange), function(x) mle[[x]]$objective)
    return(sum(mle.val))
  }
  lower <- rep(0,3)
  upper <- rep(Inf,3)
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  return(res)
}
## optimize branch lengths when force it to be ultrametric
optim.ultrametric <- function(tree,data,method="SBPLX",maxeval="100",print_level=1,...){
  ## free parameters: total tree depth, all the internal branch lengths
  ## calculate the exterior branch lengths given the free parameters
  nTips <- length(tree$tip.label)
  nIntBr <- nTips - 2 #number of interior branches
  edge <- tree$edge
  int.ind <- which(!(edge[,2] %in% c(1:nTips))) #index of interior branches in tree$edge
  ext.ind <- which(edge[,2] %in% c(1:nTips))
  fn <- function(x){ ## x <- c(intBrlen, depth)
    x <- exp(x)
    cat(x,"\n")
    depth <- x[nIntBr + 1] # tree depth is the last element in the vector
    tree$edge.length[int.ind] <- x[1:nIntBr] # assign internal branches new lengths
    tree <- ult.tree(tree,depth=depth)$tree #calculate exterior branch lengths
    cat("tree edge lengths","\n")
    cat(tree$edge.length,"\n")
    if(any(tree$edge.length<0)) return(10^7)
    else{
      res <- -mllm1(data,tree,...)$ll$loglik
      cat("loglik = ", res, "\n")
      return(res)
    }
  }
  ## constraint function: add "eval_g_ineq=g0" in the nloptr function
  #   g0 <- function(x){
  #     x <- exp(x)
  #     depth <- x[nIntBr + 1]
  #     tree$edge.length[int.ind] <- x[1:nIntBr] # assign internal branches new lengths
  #     tree <- ult.tree(tree,depth=depth) #calculate exterior branch lengths
  #     return(-tree$edge.length[ext.ind])
  #   }
  lower=rep(-Inf,nIntBr+1)
  upper=rep(Inf,nIntBr+1)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  
  res = nloptr(x0=log(c(tree$edge.length[int.ind],max(node.depth.edgelength(tree)))),
               eval_f=fn, lb=lower,ub=upper,opts=opts)
  x <- res$solution
  x <- exp(x)
  depth <- x[nIntBr + 1] # tree depth is the last element in the vector
  tree$edge.length[int.ind] <- x[1:nIntBr] # assign internal branches new lengths
  tree <- ult.tree(tree,depth=depth)$tree #calculate exterior branch lengths
  return(list(res=res,tree=tree))
}
## find internal branches of tree
# internal.br <- function(tree){
#   nTips <- length(tree$tip.label)
#   edge <- tree$edge
#   internal.ind <- which(!(edge[,2] %in% c(1:nTips)))
#   return(list(index = internal.ind,edge=edge[internal.ind,],edge.length=tree$edge.length[internal.ind]))
# }

## given a tree and desired depth, make the tree ultrametric without changing the interior branch lengths
ult.tree <- function(tree,depth=NULL){
  if(is.ultrametric(tree))
    return(list(tree=tree,depth=depth))
  else{
    nTips <- length(tree$tip.label)
    depth <- max(node.depth.edgelength(tree))
    for(i in 1:nTips){
      brs <- path_to_tip(tree,i)$br.path
      l <- length(brs)
      if(l==1)
        tree$edge.length[brs] = depth
      else
        tree$edge.length[brs[l]] = depth - sum(tree$edge.length[brs[1:(l-1)]])
    }
    return(list(tree=tree,depth=depth))
  }
}
######################################################################################################
## MLE for branch lengths, specify starting values for el, otherwise it starts with the tree supplied with branch length
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call: optim.br(data,tree,el=NULL,s=0.1,beta=be,gamma=ga,Q=NU_VEC...)
optim.br <- function(data,tree,el=NULL,method="COBYLA",maxeval="100", print_level=0, ...){
  if(is.null(el))
  {el <- tree$edge.length}
  br.num <- length(el)
  
  res.initial = mllm1(data=data,tree=tree,...)
  #these don't change with the change of s
  Qall = res.initial$Qall
  fn = function(el,data,tree){
    tree$edge.length = el
    #print(el)
    result = -mllm1(data,tree, Qall=Qall,...)$ll$loglik
    return(result)
  }
  lower=rep(0,br.num)
  upper=rep(Inf,br.num)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,
               "stopval"=-Inf,"print_level"=print_level)
  res = nloptr(x0=el,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  #print(res)
  return(res)
}

######################################################################################################

## optimize the mutation rates for nucleotides, this is only for GTR model, with 5 parameters
## can be easily changed to account for other models, e.g. Jukes-Cantor, etc.
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
optimQ <- function(tree,data,Q=rep(1,6),method="SBPLX",maxeval="100",print_level=0, ...){
  Q = Q/Q[6] #make last rate 1
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,Q=Q,...)
  #these don't change with the change of s
  fixmatall = res.initial$fixmatall
  
  ab <- Q[1:5] # optimize on the first 5 rates
  fn = function(ab,tree,data){
    #print(ab)
    result = -mllm1(data,tree,Q=c(ab,1),fixmatall=fixmatall, ...)$ll$loglik
    return(result)
  }
  lower=rep(0,5)
  upper=rep(Inf,5)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  res$solution = c(res$solution,1) # append the last rate (1) to the rate vector
  return(res)
}
######################################################################################################
## optimize parameters
## first optimize all parameters except optimal weights
## then optimize the weights of amino acids being optimal
optim.mllm1 <- function(object, optQ = FALSE, optBranch = FALSE, optsWeight = TRUE, optOpw = FALSE,optscale=FALSE,
                        control = list(epsilon=1e-08,hmaxit=10,htrace=TRUE,print_level=0,maxeval="300"),...){
  tree = object$tree
  if(any(tree$edge.length < 1e-08)){
    tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
    #object <- update(object, tree=tree)
  }
  call = object$call
  #maxit = control$maxit #maximum number of iterations for each sub optimizer
  htrace = control$htrace #print out information about steps or not?
  print_level=control$print_level
  maxeval = control$maxeval
  data = object$data
  Q = object$Q
  scale.vec = object$scale.vec
  #if(is.null(subs)) subs = c(1:(length(Q)-1),0) #default is GTR
  bfaa = object$bfaa #this is going to be the same, no matter from data or given -- empirical frequencies
  opw = NULL
  ll = object$ll$loglik
  ll1 = ll
  s = object$s
  beta = object$GMweights[2]
  gamma = object$GMweights[3]
  if(is.null(s)) s = 1
  if(is.null(beta)) beta=be
  if(is.null(gamma)) gamma=ga
  opti = TRUE # continue optimizing or not
  rounds = 0 #index of iterations
  while(opti){
    if(htrace){
      cat("\n","Round ",rounds+1,"\n")
      cat("opw = ", opw, "\n")
    }
    if(optsWeight){
      res = optim.s.weight(data,tree,maxeval=maxeval,print_level=print_level,
                           s=s,beta=beta,gamma=gamma,Q=Q,scale.vec=scale.vec,opw=opw,...)
      s = res$solution[1]
      beta = res$solution[2]
      gamma = res$solution[3]
      if(htrace){
        cat("optimize s and Grantham weights: ", ll, "--->", -res$objective, "\n")
        cat("s, beta, gamma are now: ", res$solution, "\n")
      }
      ll = -res$objective
    }
    if(optQ){
      res = optimQ(tree,data,Q=Q,maxeval=maxeval,print_level=print_level,
                   s=s,beta=beta,gamma=gamma,scale.vec=scale.vec,opw=opw, ...)
      Q = res$solution
      if(htrace){
        cat("optimize rate matrix: ", ll, "--->", -res$objective, "\n")
        cat("Q is now: ", res$solution, "\n")
      }
      ll = -res$objective
    }
    if(optBranch){
      res = optim.br(data,tree,maxeval=maxeval,print_level=print_level,
                     s=s,beta=beta,gamma=gamma,Q=Q,scale.vec=scale.vec,opw=opw, ...)
      if(htrace){
        cat("optimize branch lengths:", ll, "--->", -res$objective, "\n")
        cat("branch lengths are now: ", res$solution, "\n")
      }
      tree$edge.length = res$solution
      ll =-res$objective
    }
    if(optscale){
      res = optim.scale(data,tree,scale.vec=scale.vec,Qall=NULL,maxeval=maxeval,
                        print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q,opw=opw,...)
      scale.vec = c(1,res$solution)
      if(htrace){
        cat("optimize scales of Q:", ll, "--->", -res$objective, "\n")
        cat("scales are now: ", scale.vec, "\n")
      }
      ll = -res$objective
    }
    if(optOpw){
      res = optim.opw(data,tree,opw=opw,print_level=print_level,
                      s=s,beta=beta,gamma=gamma,Q=Q,scale.vec=scale.vec,...) #new optimizer using nloptr
      if(htrace){
        ##notice that the loglikelihood will decrease, from maximizing rule to weighted rule
        cat("optimize weights of optimal aa:", ll, "--->", -res$objective, "\n")
        cat("weights are now: ", res$solution, "\n")
      }
      opw = res$solution
      ll = -res$objective
    }
    
    rounds = rounds + 1
    if(rounds >= control$hmaxit) opti <- FALSE
    if(((ll1-ll)<=0) && ((ll1-ll)/ll < control$epsilon)) opti <- FALSE
    ll1 = ll
  }
  
  object = update(object,tree=tree,data=data,s=s,beta=beta,gamma=gamma,Q=Q,scale.vec=scale.vec,
                  opw=opw,...)
  return(object)
}
#MLE for s and Ne, using subplex method by default
#mllm1 <- function(data,tree,s=NULL,beta=be,gamma=ga,Q=NULL,dismat=NULL,fixmatall=NULL,mumat=NULL,Qall=NULL,
#                  opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,5000,beta=be,gamma=ga,Q=NU_VEC,...)
optim.s.Ne <- function(data, tree,s,Ne, method="SBPLX",maxeval="500",print_level=0,...){
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,beta=be,gamma=ga,Ne=Ne,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  
  ab <- c(s,Ne) ##initial value
  fn = function(ab,data,tree){
    cat("s, Ne = ",ab[1]," ",ab[2],"\n")
    result = -mllm1(data=data,tree=tree,s=ab[1],beta=be,gamma=ga, mumat=mumat,bfaa=bfaa,Ne=ab[2], ...)$ll$loglik
    cat(result,"\n")
    return(result)
    
  }
  lower <- rep(0,2)
  upper <- rep(Inf,2)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}

optim.Ne <- function(data, tree,s,Ne, method="SBPLX",maxeval="500",print_level=0,...){
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,beta=be,gamma=ga,Ne=Ne,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  
  ab <- Ne ##initial value
  fn = function(ab,data,tree){
    cat("Ne = ",ab,"\n")
    result = -mllm1(data=data,tree=tree,s=s,beta=be,gamma=ga, mumat=mumat,bfaa=bfaa,Ne=ab, ...)$ll$loglik
    cat(result,"\n")
    return(result)
    
  }
  lower <- 0
  upper <- Inf
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}