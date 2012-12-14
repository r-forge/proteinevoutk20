subplex2.opw <- function(data, tree,opw=NULL, bad.value=-10000000, ...){
  if(is.null(opw))
    opw = findBf2(data)
  opw = opw/sum(opw)
  l = length(opw)
#   opw[opw==0] <- 1e-8
#   lopw = log(opw)
  #nenner = 1/opw[l]
  #lopw = log(opw*nenner) #scale the vector by the last entry
   # optimize on the all entries except the last one
  
  res = mllm(data=data,tree=tree,opw=opw,...) #store llmat (loglikelihood values for all opaa) 
  llmat = exp(res$ll$llmat)    #so that they don't need to be evaluated again and again
  weight = attr(data,"weight")
  print(res$ll$loglik)
  
  #opw = opw[-l]
  
#   fn = function(opw){
#     #opw <- append(opw, 1-sum(opw))
#     #opw=opw/sum(opw)
#     result<-NA
#     if(min(opw)<0 || max(opw)>1) {
#       result <- bad.value
#     }
#     else {
#       sitelik = log(llmat %*% opw)
#       result = sum(sitelik*weight)
#     }
#     cat("par:",opw,"val:",result,"\n")
#     
#     return(-result)
#   }
  eval_f_list <- function(opw){
    return(list("objective"=llaaw1(opw,weight=weight,llmat=llmat),
                "gradient"=llaaw_grad(opw,weight=weight,llmat=llmat)))
  }
  eval_g_list <- function(opw){
    return(list("constraints"=sum(opw)-1,"jacobian"=rep(1,20)))
  }
  lower <- rep(0,20)
  upper <- rep(1,20)
  #ui=rep(-1,19)
  #ci=-1
  #res <- constrOptim(theta=opw,f=fn,grad=f_grad,ui=ui,ci=ci)
  local_opts <- list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1e-7)
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","maxeval"="1000000","xtol_rel"=1e-7,"ftol_rel"=.Machine$double.eps,"local_opts"=local_opts)
  res = nloptr(x0=opw,eval_f=eval_f_list,eval_g_eq=eval_g_list, lb=lower,ub=upper,opts=opts)
  #print(res[[2]])
#   opw = exp(res$solution)
#   opw <- append(opw,1-sum(opw))
#   res$par = opw
  return(res)
}

subplex.opw <- function(data, tree,opw=NULL, ...){
  if(is.null(opw))
    opw = findBf2(data)
  l = length(opw)
  nenner = 1/opw[l]
  lopw = log(opw*nenner) #scale the vector by the last entry
  lopw = lopw[-l] # optimize on the all entries except the last one
  res = mllm(data=data,tree=tree,opw=opw,...)
  exp_mat = exp(res$ll$llmat)
  weight = attr(data,"weight")
  
  fn = function(lopw){
    opw = exp(c(lopw,0))
    opw=opw/sum(opw)
    sitelik = log(exp_mat %*% opw)
    result = sum(sitelik*weight)
    cat("par:",opw,"val:",result,"\n")
    return(-result)
  }
  lower <- rep(-1000,19)
  upper <- rep(100,19)
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="1000000","ftol_rel"=.Machine$double.eps)
  res = nloptr(x0=lopw,eval_f=fn, lb=lower,ub=upper,opts=opts)
  #print(res[[2]])
  opw = exp(c(res[[1]],0))
  opw = opw/sum(opw)
  res$par = opw
  return(res)
}
