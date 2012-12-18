optim.opw <- function(data, tree,opw=NULL, bad.value=-10000000, ...){
  if(is.null(opw))
    opw = findBf2(data)
  opw = opw/sum(opw) #opw given in the function call doesn't have to sum to 1
  
  res = mllm(data=data,tree=tree,opw=opw,...) #store llmat (loglikelihood values for all opaa) 
  llmat = exp(res$ll$llmat)    #so that they don't need to be evaluated again and again
  weight = attr(data,"weight")
  print(res$ll$loglik) #function value at the starting point
  
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
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","maxeval"="1000000","xtol_rel"=1e-7,"ftol_rel"=.Machine$double.eps,"local_opts"=local_opts)
  res = nloptr(x0=opw,eval_f=eval_f_list,eval_g_eq=eval_g_list, lb=lower,ub=upper,opts=opts)
  # res$objective: best function value found
  # res$solution: best parameter values
  return(res)
}