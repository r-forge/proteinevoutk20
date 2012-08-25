source("~/proteinevoutk20/pkg/R/simulation.R")

find_hessian <- function(mle,tree,data,alpha,MuMat,m=20,protein_op=NULL,root=NULL,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.37e07){
  negloglike <- function(para){
    s = para[1]
    beta = para[2]
    gamma = para[3]
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,C,Phi,q,Ne))
}
  return(hessian(negloglike,mle,method="Richardson"))
}

find_hessian_log<- function(mle,tree,data,alpha,MuMat,m=20,protein_op=NULL,root=NULL,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.37e07){
  negloglike <- function(para){
    print(para)
    s = exp(para[1])
    beta = exp(para[2])
    gamma = exp(para[3])
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,C,Phi,q,Ne))
}
  return(hessian(negloglike,mle,method="Richardson"))
}
