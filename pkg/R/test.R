## MLE for opw (weights for optimal amino acids)
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call: optim.opw.subplex(data,tree,s=0.1,Q=NU_VEC...)
optim.opw.subplex <- function(data, tree,opw=rep(1/20,20), ...){
  l = length(opw)
  nenner = 1/opw[l]
  lopw = log(opw*nenner) #scale the vector by the last entry
  lopw = lopw[-l] # optimize on the all entries except the last one
  fn = function(lopw,data,tree, ...){
    opw = exp(c(lopw,0))
    opw=opw/sum(opw)
    result = mllm(data=data,tree=tree,opw=opw, ...)$ll$loglik
    return(result)
  }
  res = subplex(par=lopw,fn=fn,data=data,tree=tree,...)
  #res = optim(par=lopw,fn=fn,gr=NULL,method=method,lower=-Inf,upper=Inf,
              #control=list(fnscale=-1,trace=trace,maxit=maxit),data=data,tree=tree, ...)
  #print(res[[2]])
  opw = exp(c(res[[1]],0))
  opw = opw/sum(opw)
  res$par = opw
  return(res)
}