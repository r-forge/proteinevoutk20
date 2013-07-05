###########################################################################
## look at the nonzero entries in the AArate matrix, and try to figure out the 
## influence of NU_rates on AA_rates
# nonzeroPair <- indAApair[!sapply(indAApair,FUN=is.null)]
# sapply(1:length(nonzeroPair),function(x) nonzeroPair[[x]]$rates)
# sapply(1:length(nonzeroPair),function(x) nonzeroPair[[x]]$tot)
# NU_to_AA <- matrix(0,nrow=6,ncol=length(nonzeroPair))
# for(i in 1:length(nonzeroPair)){
#   rates <- nonzeroPair[[i]]$rates
#   NU_to_AA[,i] <- sapply(1:6,function(x) sum(rates==x))/nonzeroPair[[i]]$tot
# }
# get.root <- function(mat){
#   l <- dim(mat)[2]
#   sapply(1:l,function(x) which(mat[,x]==1))
# }
optimQ <- function(tree,data,Q=rep(1,6),method="SBPLX",maxeval="100",print_level=0, ...){
  Q = Q/Q[6] #make last rate 1
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,Q=Q,...)
  #these don't change with the change of s
  fixmatall = res.initial$fixmatall
  
  ab <- Q[1:5] # optimize on the first 5 rates
  fn = function(ab,tree,data){
    cat(ab,"\n")
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

hessianQ <- function(mle,data,tree,...){
  fn <- function(Q){
    print(Q)
    return(-mllm1(data=data,tree=tree,Q=Q,...)$ll$loglik)
  }
  return(hessian(fn,mle,method="Richardson"))
}