#source("~/proteinevoutk20/pkg/R/simulation.R")
# find the hessian matrix at the mle for c(s,beta,gamma) for the -loglikelihood function
# since it's the -loglikelihood, the covariance matrix will just be the inverse of hessian
# instead of the inverse of the -hessian
# square roots of diagonals from the inverse hessian will give the standard errors
source("~/proteinevoutk20/pkg/R/rokas_collect.R")
source("~/proteinevoutk20/pkg/R/main.R")
find_hessian <- function(mle,data,tree,...){
  fn <- function(para){
    print(para)
    s = para[1]
    beta = para[2]
    gamma = para[3]
    return(-mllm1(data,tree,s=s,beta=beta,gamma=gamma,...)$ll$loglik)
  }
  return(hessian(fn,mle,method="Richardson"))
}
## After sourcing the file "rokas_collect" to get the list res_max
hes <- vector("list",length=106)
sd_max <- matrix(0,nrow=3,ncol=106)
for(i in 1:106){
  res <- res_max[[i]]
  hes_mat <- find_hessian(c(res$s,res$GMweights[2:3]),res$data,res$tree,Q=res$Q,bfaa=res$bfaa,opaa=res$ll$opaa)
  std <- sqrt(diag(solve(hes_mat)))
  hes[[i]] <- hes_mat
  sd_max[,i] <- std
}
save(hes,sd_max,res_max,file="rokas_hessian.RData")