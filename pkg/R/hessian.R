#source("~/proteinevoutk20/pkg/R/simulation.R")
# find the hessian matrix at the mle for c(s,beta,gamma) for the -loglikelihood function
find_hessian <- function(mle,data,tree,...){
  fn <- function(para){
    print(para)
    s = para[1]
    beta = para[2]
    gamma = para[3]
    return(-mllm(data,tree,s=s,beta=beta,gamma=gamma,...)$ll$loglik)
  }
  return(hessian(fn,mle,method="Richardson"))
}
