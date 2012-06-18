#substitution rates calculated according to S-H formula
rate_mat <- function(fit_vec, Ne, mu){
  d <- length(fit_vec) # number of states in the system
  mat <- matrix(NA,nrow=d,ncol=d)
  for(i in 1:d){
    for(j in 1:d){
      if(j!=i){
        if(fit_vec[i]!=fit_vec[j]){
          pi_ij <- (1-fit_vec[i]/fit_vec[j])/(1-(fit_vec[i]/fit_vec[j])^(2*Ne))
          mat[i,j] <- 2*Ne*mu*pi_ij
        }
        else{
          mat[i,j] <- mu
        }
      } #end if
    }#end for j
    mat[i,i] <- -sum(mat[i,],na.rm=TRUE)
  }# end for i
  return(mat)
} #end function

#substitution calculated according to canonical formula
rate_mat_cnl <- function(fit_vec, Ne, mu){
  d <- length(fit_vec) # number of states in the system
  mat <- matrix(NA,nrow=d,ncol=d)
  for(i in 1:d){
    for(j in 1:d){
      if(j!=i){
        if(fit_vec[i]!=fit_vec[j]){
          s <- (fit_vec[j]-fit_vec[i])/fit_vec[i] #selection advantage of j comparing to i
              #(w_j - w_i)/w_i
          pi_ij <- (1-exp(-s))/(1-exp(-2*Ne*s))
          mat[i,j] <- 2*Ne*mu*pi_ij
        }
        else{
          mat[i,j] <- mu
        }
      } #end if
    }#end for j
    mat[i,i] <- -sum(mat[i,],na.rm=TRUE)
  }# end for i
  return(mat)
} #end function
fit_vec <- c(9.99999,9.999991,9.999992,9.99999) #fitness vector
Ne <- 1000 #effective population size##
mu <- 0.1 #mutation rates are all the same for now
W <- rate_mat(fit_vec,Ne,mu)
W_cnl <- rate_mat_cnl(fit_vec,Ne,mu)
freq_t <- expm(W*10000)[1,] #freq <- fit_vec^(2*Ne-1)/sum(fit_vec^(2*Ne-1))
freq <- (fit_vec/max(fit_vec))^(2*Ne-1) #this is a better way to evaluate S-H formula
freq <- freq/sum(freq)
sprintf("%.20f",freq-freq_t)

#This is to verify that the two formulae are close enough when s << 1
s <- seq(0,0.1,0.0001)
pi_1 <- (1-exp(-s))/(1-exp(-2*Ne*s))
pi_2 <- (1-(1/(1+s)))/(1-(1/(1+s))^(2*Ne))
plot(pi_1,type='l')
points(pi_2,col="red",type="l")
