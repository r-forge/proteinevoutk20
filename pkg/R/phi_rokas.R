# 
# 
# lik_phi <- function(phi_ij,phi_hat,sigma,sigma_hat){ 
#   #switch sigma and sigma_hat doesn't change the result
#   integrand <- function(phi_i)
#     #(1/(2*pi*sigma*sigma_hat))*exp(-(phi_ij-phi_i)^2/(2*sigma^2)-(phi_i-phi_hat)^2/(2*sigma_hat^2))
#     dnorm(phi_ij,mean=phi_i,sd=sigma)*dnorm(phi_i,mean=phi_hat,sd=sigma_hat)
#   integrate(integrand,-Inf,Inf)$value
# }
# 
# lik_phi_vec <- Vectorize(lik_phi,vectorize.args="phi_ij")
# optim.phi.sigma <- function(phis,phi_hat,sigma,sigma_hat,maxeval="500",print_level=1, ...){
#   ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
#   ab <- c(phi_hat,sigma,sigma_hat)
#   fn <- function(ab,phis){
#     cat(ab,"\n")
#     loglik <- -sum(log(lik_phi_vec(phis,ab[1],ab[2],ab[3])))
#     return(loglik) #summation of all values
#   }
#   
#   lower <- c(-Inf,0.01,0.01) #lower bound
#   upper <- c(Inf,Inf,Inf) #upper bound
#   #options for optimizer
#   opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"=maxeval,
#                "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
#   res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,phis=phis)
#   return(res)
# }
# #likelihood function of log(g*phi), 
# #given mean, variance value of phi, a, b, and variance of error sigma_e
# lik_gphi <- function(gphi,phi,a,b,sigma,sigma_e){
#   dnorm(gphi,mean=a+(b+1)*phi,sd=sqrt(((b+1)^2*sigma^2)+sigma_e^2))
# }
# lik_gphi_vec <- Vectorize(lik_gphi,vectorize.args=c("gphi","phi"))
# 
# lik_gphi_fix <- function(gphi,phi,a,b,c){
#   dnorm(gphi,mean=a+(b+1)*phi,sd=c)
# }
# lik_gphi_fix_vec <- Vectorize(lik_gphi_fix,vectorize.args=c("gphi","phi"))
# 
# optim.ab.sigma_e <- function(gphis,phis, a,b,c,maxeval="500",print_level=1){
#   ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
#   ab <- c(a,b,log(c))
#   fn <- function(ab){
#     cat(ab[1],ab[2],exp(ab[3]),"\n")
#     loglik <- -sum(log(lik_gphi_fix_vec(gphi=gphis,phi=phis,a=ab[1],b=ab[2],
#                                     c=exp(ab[3]))))
#     return(loglik) #summation of all values
#   }
#   
#   lower <- c(-Inf,-Inf,-Inf) #lower bound
#   upper <- c(Inf,Inf,Inf) #upper bound
#   #options for optimizer
#   opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"=maxeval,
#                "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
#   res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
#   return(res)
# }
## phi is a vector, represent observed values of phi for ith gene
library(nloptr)
source("phi.R")
lik_y_phi <- function(gphi,phi,mu,a,b,phi.sd,mu.sd,phig.sd){
  integrand <- function(phi_i)
    prod(dnorm(phi_i,mu,mu.sd)*dnorm(gphi,a+(b+1)*phi_i,phig.sd)*dnorm(phi,phi_i,phi.sd))
  integrand_vec <- Vectorize(integrand,vectorize.args="phi_i")
  integrate(integrand_vec,-Inf,Inf)$value
}
lik_y_phi_vec <- function(gphi_vec,phi_mat,mu,a,b,phi.sd,mu.sd,phig.sd){
  phi_list <- lapply(1:dim(phi_mat)[1], function(i) as.numeric(phi_mat[i,]))
  mapply(lik_y_phi,gphi_vec,phi_list,
         MoreArgs=list(mu=mu,a=a,b=b,phi.sd=phi.sd,mu.sd=mu.sd,phig.sd=phig.sd))

}
optim.ab.sd <- function(gphis,phis, mu,a,b,phi.sd,mu.sd,phig.sd,maxeval="500",print_level=1){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  ab <- c(mu,a,b,log(phi.sd),log(mu.sd),log(phig.sd))
  fn <- function(ab){
    cat(ab[1],ab[2],ab[3],exp(ab[4]),exp(ab[5]),exp(ab[6]),"\n")
    loglik <- -sum(log(lik_y_phi_vec(gphi_vec=gphis,phi_mat=phis,mu=ab[1],a=ab[2],b=ab[3],
                                        phi.sd=exp(ab[4]),mu.sd=exp(ab[5]),phig.sd=exp(ab[6]))))
    return(loglik) #summation of all values
  }
  
  lower <- rep(-100,6) #lower bound
  upper <- rep(100,6) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"=maxeval,
               "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  return(res)
}

est.logphi <- data.matrix(phidata[,1:3])
est.logphig <- lgphi
mu=-1
a=1
b=1
phi.sd=1
mu.sd=1
phig.sd=1
print(-sum(log(lik_y_phi_vec(est.logphig,est.logphi,mu,a,b,phi.sd,mu.sd,phig.sd))))
opt.para <- optim.ab.sd(est.logphig,est.logphi,-2.323789, 1.825016, -0.5452444, exp(0.673978), exp(-4.971654), exp(-1.007252))
# cov.wy <- cov(est.logphi,est.logphig)
# cov.wy1 <- cov(est.logphi1,est.logphig)
# var.w <- var(est.logphi)
# var.w1 <- var(est.logphi1)
# s.u <- var(est.logphi-est.logphi1)/2
# b1 <- cov.wy/(var.w1-s.u)
# b2 <- cov.wy1/(var.w1-s.u)
# b3 <- cov.wy1/(var.w-s.u)
