## MLE for opw (weights for optimal amino acids)
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call: optim.opw.subplex(data,tree,s=0.1,Q=NU_VEC...)
optim.opw.new <- function(data, tree,opw=rep(1/20,20),method="Nelder-Mead", maxit=500, trace=0, ...){
  l = length(opw)
  nenner = 1/opw[l]
  lopw = log(opw*nenner) #scale the vector by the last entry
  lopw = lopw[-l] # optimize on the all entries except the last one
  res = mllm(data,tree,...)
  llmat = res$ll$llmat
  weight = attr(res$data,"weight")
  fn = function(lopw,data,tree, ...){
    opw = exp(c(lopw,0))
    opw=opw/sum(opw)
    result = weight %*% llmat %*% opw
    cat("par:",opw,"val:",result,"\n")
    return(result)
  }
  res = optim(par=lopw,fn=fn,gr=NULL,method=method,lower=-Inf,upper=Inf,
              control=list(fnscale=-1,trace=trace,maxit=maxit),data=data,tree=tree, ...)
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

source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/readRokas.R")
tree <- ROKAS_TREE
gene2 <- ROKAS_DATA[[2]]
#starting point of optimization
opw.s <- rep(0.001,20)
opw.s[1] <- 1
opw.s
opw = optim.opw(gene2,tree,opw=opw.s,maxit=2500,trace=1,s=0.1,Q=NU_VEC)
opw

#loglikelihoods for gene2 when opaa is a particular amino acid (20 by 1 vector)
ll <- NULL
for(i in 1:20){
  opw = rep(0,20)
  opw[i] = 1
  cat(i,":",opw, "\n")
  res = mllm(ROKAS_DATA[[2]],tree,0.1,be,ga,Q=NU_VEC,opw=opw)
  ll = c(ll,res$ll$loglik)
}

#collect results on optimization of opw for gene2, gene1 can be attained similarly
dir <- "~/proteinevoutk20/pkg/scratch/lab9/RokasOpw/Gene2/"
res.opw <- vector("list",length=20)
l <- 20
for(genect in 1:l){
  filename = paste(dir,"gene2_",genect,".RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res.opw[[genect]] <- opw
}

opw_par <- sapply(1:20,function(x) res.opw[[x]]$par)
opw_value <- sapply(1:20,function(x) res.opw[[x]]$value)
opw_count <- sapply(1:20,function(x) res.opw[[x]]$counts[1])

#for a given vector of log likelihood, and the corresponding opaa, find the smallest set of aa
# that cover the 95% of the total likelihood
aa.set <- function(llvec){
  lvec = exp(llvec) #likelihood
  ord = order(llvec,decreasing=T) #order of llvec (increasing)
  lvec = lvec[ord] # ordered likelihood
  tol = sum(lvec) #total likelihood
  CumSum = cumsum(lvec) # cumulative sum 
  ind = sum(CumSum < tol*0.95) + 1 #the index that gives > 95% totol likelihood
  Sum = CumSum[ind] # sum of the first "ind" terms
  #return the set of amino acids that give > 95% likelihood, percentile, and 
  #the percentile of the likelihood given by the optimal aa (max rule)
  return(list(aa=ord[1:ind], percentile=Sum/tol,op.percentile=lvec[1]/tol))
}