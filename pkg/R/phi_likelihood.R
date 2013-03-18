lik_phi <- function(phi_ij,phi_hat,sigma,sigma_hat){ 
  #switch sigma and sigma_hat doesn't change the result
  integrand <- function(phi_i)
    #(1/(2*pi*sigma*sigma_hat))*exp(-(phi_ij-phi_i)^2/(2*sigma^2)-(phi_i-phi_hat)^2/(2*sigma_hat^2))
    dnorm(phi_ij,mean=phi_i,sd=sigma)*dnorm(phi_i,mean=phi_hat,sd=sigma_hat)
  integrate(integrand,-Inf,Inf)$value
}

lik_phi_vec <- Vectorize(lik_phi,vectorize.args="phi_ij")
optim.phi.sigma <- function(phis,phi_hat,sigma,sigma_hat,maxeval="500",print_level=1, ...){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  ab <- c(phi_hat,sigma,sigma_hat)
  fn <- function(ab,phis){
    cat(ab,"\n")
    loglik <- -sum(log(lik_phi_vec(phis,ab[1],ab[2],ab[3])))
    return(loglik) #summation of all values
  }
  
  lower <- c(-Inf,0.01,0.01) #lower bound
  upper <- c(Inf,Inf,Inf) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"=maxeval,
               "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,phis=phis)
  return(res)
}

lik_gphi <- function(gphi,phi,sigma,a,b,sigma_e){
  dnorm(gphi,mean=a+(b+1)*phi,sd=sqrt(((b+1)^2*sigma^2)+sigma_e^2))
}
lik_gphi_vec <- Vectorize(lik_gphi,vectorize.args=c("gphi","phi"))
optim.ab.sigma_e <- function(gphis,phis, a,b,sigma,sigma_e,maxeval="500",print_level=1){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  ab <- c(a,b,log(sigma_e))
  fn <- function(ab){
    cat(ab,"\n")
    loglik <- -sum(log(lik_gphi_vec(gphi=gphis,phi=phis,sigma=sigma,a=ab[1],b=ab[2],
                                    sigma_e=exp(ab[3]))))
    return(loglik) #summation of all values
  }
  
  lower <- c(-100,-100,-Inf) #lower bound
  upper <- c(100,100,Inf) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"=maxeval,
               "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  return(res)
}

#ps1 <- optim.phi.sigma(allphi,1,1,1)
#absigma <- optim.ab.sigma_e(lgphi,a=0,b=0,phi_hat=ps1$solution[1],sigma_hat=ps1$solution[3],sigma_e=1)
