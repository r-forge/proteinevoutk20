optr.w <- function(beta,gamma,generange,tree,multicore=FALSE,...){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  ab <- c(beta,gamma)
  fn <- function(ab,generange,tree){
    #print(ab)
    mle <- optim.s.range(ab[1],ab[2],generange,tree=tree,multicore=multicore,...) #call the previous function to optimize s for all genes
    mle.val <- sapply(1:length(generange),function(ind) mle[[ind]]$objective) # best -loglikelihood values
    return(sum(mle.val)) #summation of all values
  }
  
  lower <- c(0,0) #lower bound
  upper <- c(Inf,Inf) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-7,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,generange=generange,tree=tree)
  return(res)
}
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,trace=1,Q=NU_VEC))
optr.s.weight <- function(data, tree, s,beta,gamma, ...){
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,beta=beta,gamma=gamma,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  
  ab <- c(s,beta,gamma) ##initial value
  fn = function(ab,data,tree){
    #print(ab)
    result = -mllm1(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], mumat=mumat,bfaa=bfaa, ...)$ll$loglik
    return(result)
  }
  lower <- rep(0,3)
  upper <- rep(Inf,3)
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-5,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}

optimQ <- function(tree,data,Q=rep(1,6), ...){
  Q = Q/Q[6] #make last rate 1
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,Q=Q,...)
  #these don't change with the change of s
  fixmatall = res.initial$fixmatall
  bfaa = res.initial$bfaa
  
  ab <- Q[1:5] # optimize on the first 5 rates
  fn = function(ab,tree,data){
    print(ab)
    result = -mllm1(data,tree,Q=c(ab,1),fixmatall=fixmatall,bfaa=bfaa, ...)$ll$loglik
    return(result)
  }
  lower=rep(0,5)
  upper=rep(Inf,5)
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-5,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  res$solution = c(res$solution,1) # append the last rate (1) to the rate vector
  #print(res)
  return(res)
}

optimBr <- function(data,tree,el=NULL, ...){
  if(is.null(attr(tree,"order")) || attr(tree,"order") == "cladwise")
    tree <- reorderPruning(tree)
  if(is.null(el))
    {el <- tree$edge.length}
  br.num <- length(el)
  
  res.initial = mllm1(data=data,tree=tree,...)
  #these don't change with the change of s
  Qall = res.initial$Qall
  bfaa = res.initial$bfaa
  fn = function(el,data,tree){
    tree$edge.length = el
    print(el)
    result = -mllm1(data,tree, Qall=Qall,bfaa=bfaa,...)$ll$loglik
    return(result)
  }
  lower=rep(0,br.num)
  upper=rep(Inf,br.num)
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-5,"ftol_rel"=.Machine$double.eps^0.5,
               "stopval"=-Inf,"print_level"=1)
  res = nloptr(x0=el,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  print(res)
  return(res)
}