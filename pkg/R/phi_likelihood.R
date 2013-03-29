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
dlik <- function(gphi,phi,phi_i,mu,a,b,phi.sd,mu.sd,phig.sd)
  prod(dnorm(phi_i,mu,mu.sd)*dnorm(gphi,a+(b+1)*phi_i,phig.sd)*dnorm(phi,phi_i,phi.sd))

dlik_log <- function(gphi,phi,phi_i,mu,a,b,phi.sd,mu.sd,phig.sd){
  dlik.log <- sum(dnorm(phi_i,mu,mu.sd,log=T)+dnorm(gphi,a+(b+1)*phi_i,phig.sd,log=T)+dnorm(phi,phi_i,phi.sd,log=T))
  return(dlik.log)
}
lik_y_phi <- function(gphi,phi,mu,a,b,phi.sd,mu.sd,phig.sd){
  integrand <- function(phi_i)
    prod(dnorm(phi_i,mu,mu.sd)*dnorm(gphi,a+(b+1)*phi_i,phig.sd)*dnorm(phi,phi_i,phi.sd))
  integrand_vec <- Vectorize(integrand,vectorize.args="phi_i")
  integrate(integrand_vec,-Inf,Inf)$value
}
lik_y_phi_vec <- function(gphi_vec,phi_mat,mu,a,b,phi.sd,mu.sd,phig.sd){
  phi_list <- lapply(1:dim(phi_mat)[1], function(i) as.numeric(phi_mat[i,])) #convert matrix to list
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
  
  lower <- rep(-10,6) #lower bound
  upper <- rep(10,6) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=maxeval,
               "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  return(res)
}
#ps1 <- optim.phi.sigma(allphi,1,1,1)
#absigma <- optim.ab.sigma_e(lgphi,a=0,b=0,phi_hat=ps1$solution[1],sigma_hat=ps1$solution[3],sigma_e=1)
#nreps<-20
ngenes<-106
#true.vals<-matrix(ncol=3, nrow=0)
#est.vals<-matrix(ncol=3, nrow=0)
phi.sd <- 1
phig.sd <- 0.01
mu.sd <- 0.5
a<-3
b<-0
#c<-runif(1, 0.000, 0.01)
mu <- runif(1,-3,3)
#true.vals <- rbind(true.vals, c(a, b, c))
#true.logphi <- runif(ngenes, -6, -1)
true.logphi <- rnorm(ngenes, mu,mu.sd)
est.logphi1 <- rnorm(length(true.logphi), true.logphi, phi.sd)
est.logphi2 <- rnorm(length(true.logphi), true.logphi, phi.sd)
est.logphi3 <- rnorm(length(true.logphi), true.logphi, phi.sd)
est.logphi <- cbind(est.logphi1,est.logphi2,est.logphi3)
est.logphig <- rnorm(length(true.logphi), a+(b+1)*true.logphi,phig.sd)
par(mfrow=c(2,3))
plot(true.logphi,est.logphi1)
plot(est.logphi2~true.logphi)
plot(est.logphi3~true.logphi)
plot((est.logphig-true.logphi)~true.logphi)
plot((est.logphig-est.logphi1)~est.logphi1)
plot((est.logphig-est.logphi2)~est.logphi2)
lm <- lm((est.logphig-true.logphi)~true.logphi)
lm1 <- lm((est.logphig-est.logphi1)~est.logphi1)
lm2 <- lm((est.logphig-est.logphi2)~est.logphi2)

#plot(est.logphig-est.logphi2~est.logphi2)
#print(-sum(log(lik_y_phi_vec(est.logphig,est.logphi,mu,a,b,phi.sd,mu.sd,phig.sd))))
#optim.ab.sd(est.logphig,est.logphi,mu,a,b,phi.sd,mu.sd,phig.sd)
# cov.wy <- cov(est.logphi,est.logphig)
# cov.wy1 <- cov(est.logphi1,est.logphig)
# var.w <- var(est.logphi)
# var.w1 <- var(est.logphi1)
# s.u <- var(est.logphi-est.logphi1)/2
# b1 <- cov.wy/(var.w1-s.u)
# b2 <- cov.wy1/(var.w1-s.u)
# b3 <- cov.wy1/(var.w-s.u)
