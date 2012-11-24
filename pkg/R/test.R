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
    cat(opw,"\n")
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

#MLE for s, beta and gamma, using Nelder-Mead method by default
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,Q=NU_VEC)
optim.s.weight.pso <- function(data, tree, s,beta,gamma, Q){
  ab <- c(s,beta,gamma)
#   ab[ab==0] <- 1e-08 #take care of log(0)
#   ab <- log(ab)
  fn = function(ab){
#     ab = exp(ab)
    result = mllm(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], Q=Q)$ll$loglik
    cat("par:",ab,"val:",result,"\n")
    return(result)
  }
  res = optim_pso(objective_function=fn,number_of_parameters=3,number_of_particles=4,max_number_of_iterations=100,
                  parameter_bounds=cbind(rep(0,3),rep(10,3)),initial_estimates=matrix(ab,ncol=1),logfile="ppso.log",load_projectfile="no",
                  do_plot="base")
  res$par = exp(res$par)
  return(res)
}

f <- function(sd, Ne){

  if(sd==0)
    return(1/(2*Ne))
  else{
    sd = mpfr(sd,prec=10000)
    Ne = mpfr(Ne,prec=10000)
    return((1-exp(sd))/(1-exp(sd*2*Ne)))
  }
}