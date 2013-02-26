#source("~/proteinevoutk20/pkg/R/simulation.R")
# find the hessian matrix at the mle for c(s,beta,gamma) for the -loglikelihood function
# since it's the -loglikelihood, the covariance matrix will just be the inverse of hessian
# instead of the inverse of the -hessian
# square roots of diagonals from the inverse hessian will give the standard errors
find_hessian <- function(mle,data,tree,...){
  fn <- function(para){
    #print(para)
    s = para[1]
    beta = para[2]
    gamma = para[3]
    return(-mllm1(data,tree,s=s,beta=beta,gamma=gamma,...)$ll$loglik)
  }
  return(hessian(fn,mle,method="Richardson"))
}
