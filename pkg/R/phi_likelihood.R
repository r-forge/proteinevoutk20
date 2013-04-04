# library(nloptr)
# dlik <- function(gphi,phi,phi_i,mu,a,b,phi.sd,mu.sd,phig.sd)
#   prod(dnorm(phi_i,mu,mu.sd)*dnorm(gphi,a+(b+1)*phi_i,phig.sd)*dnorm(phi,phi_i,phi.sd))
# 
# dlik_log <- function(gphi,phi,phi_i,mu,a,b,phi.sd,mu.sd,phig.sd){
#   dlik.log <- sum(dnorm(phi_i,mu,mu.sd,log=T)+dnorm(gphi,a+(b+1)*phi_i,phig.sd,log=T)+dnorm(phi,phi_i,phi.sd,log=T))
#   return(dlik.log)
# }
# lik_y_phi <- function(gphi,phi,mu,a,b,phi.sd,mu.sd,phig.sd){
#   integrand <- function(phi_i)
#     prod(dnorm(phi_i,mu,mu.sd)*dnorm(gphi,a+(b+1)*phi_i,phig.sd)*dnorm(phi,phi_i,phi.sd))
#   integrand_vec <- Vectorize(integrand,vectorize.args="phi_i")
#   integrate(integrand_vec,-Inf,Inf)$value
# }
# lik_y_phi_vec <- function(gphi_vec,phi_mat,mu,a,b,phi.sd,mu.sd,phig.sd){
#   phi_list <- lapply(1:dim(phi_mat)[1], function(i) as.numeric(phi_mat[i,])) #convert matrix to list
#   mapply(lik_y_phi,gphi_vec,phi_list,
#          MoreArgs=list(mu=mu,a=a,b=b,phi.sd=phi.sd,mu.sd=mu.sd,phig.sd=phig.sd))
# 
# }
# optim.ab.sd <- function(gphis,phis, mu,a,b,phi.sd,mu.sd,phig.sd,maxeval="500",print_level=1){
#   ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
#   ab <- c(mu,a,b,phi.sd,mu.sd,phig.sd)
#   fn <- function(ab){
#     cat(ab[1],ab[2],ab[3],ab[4],ab[5],ab[6],"\n")
#     loglik <- -sum(log(lik_y_phi_vec(gphi_vec=gphis,phi_mat=phis,mu=ab[1],a=ab[2],b=ab[3],
#                                         phi.sd=ab[4],mu.sd=ab[5],phig.sd=ab[6])))
#     return(loglik) #summation of all values
#   }
#   
#   lower <- c(-10,-10,-10,0,0,0) #lower bound
#   upper <- c(10,10,10,10,10,10) #upper bound
#   #options for optimizer
#   opts <- list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=maxeval,
#                "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
#   res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
#   return(res)
# }
# 
# optim.ab <- function(gphis,phis, mu,a,b,phi.sd,mu.sd,phig.sd,maxeval="500",print_level=1){
#   ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
#   ab <- c(mu,a,b)
#   fn <- function(ab){
#     cat(ab[1],ab[2],ab[3],"\n")
#     loglik <- -sum(log(lik_y_phi_vec(gphi_vec=gphis,phi_mat=phis,mu=ab[1],a=ab[2],b=ab[3],
#                                      phi.sd=phi.sd,mu.sd=mu.sd,phig.sd=phig.sd)))
#     return(loglik) #summation of all values
#   }
#   
#   lower <- c(-10,-10,-10) #lower bound
#   upper <- c(10,10,10) #upper bound
#   #options for optimizer
#   opts <- list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=maxeval,
#                "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
#   res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
#   return(res)
# }
#ps1 <- optim.phi.sigma(allphi,1,1,1)
#absigma <- optim.ab.sigma_e(lgphi,a=0,b=0,phi_hat=ps1$solution[1],sigma_hat=ps1$solution[3],sigma_e=1)
#nreps<-20
ngenes<-106
#true.vals<-matrix(ncol=3, nrow=0)
#est.vals<-matrix(ncol=3, nrow=0)
phi.sd <- 2
phig.sd <- 0.1
mu.sd <- 0.5
a<-3
b<-2
#c<-runif(1, 0.000, 0.01)
mu <- 10
#true.vals <- rbind(true.vals, c(a, b, c))
#true.logphi <- runif(ngenes, -6, -1)
true.logphi <- rnorm(ngenes, mu,mu.sd)
#true.logphi <- log(rokasPhi$SEMPPR)
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
gphi.lm <- lmodel2(est.logphig~est.logphi1,range.y="interval",range.x="interval",nperm=99)
true.gphi.lm <- lmodel2(est.logphig~true.logphi,range.y="interval",range.x="interval",nperm=99)

# plot((est.logphig)~true.logphi)
# plot((est.logphig)~est.logphi1)
# plot((est.logphig)~est.logphi2)

