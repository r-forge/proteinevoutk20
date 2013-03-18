

phi.sd <- 0.05
phig.sd <- 0.01


GetLikelihoodSingle <- function(logphi, logphig, a, b, c) {
   return(dnorm(logphig, mean=a+(b+1)*logphi, sd=max(0,sqrt(c), na.rm=TRUE), log=TRUE))
}

GetLikelihoodMany <- function(logphi.vector, logphig.vector, a, b, c) {
   neglnL<-0
   for (i in sequence(length(logphi.vector))) {
      neglnL <- neglnL - GetLikelihoodSingle(logphi.vector[i], logphig.vector[i], a, b, c) 
   }
   return(neglnL)
}

LikelihoodForOptim <- function(x, logphi.vector, logphig.vector) {
  return(GetLikelihoodMany(logphi.vector=logphi.vector, logphig.vector=logphig.vector, a=x[1], b=x[2], c=exp(x[3]))) 
}

nreps<-2
ngenes<-1060
true.vals<-matrix(ncol=3, nrow=0)
est.vals<-matrix(ncol=3, nrow=0)
for (rep in sequence(nreps)) {
  a<-runif(1, -1, 3)
  b<-runif(1, -3, 3)
  c<-runif(1, 0.000, 0.01)
  true.vals <- rbind(true.vals, c(a, b, c)) #store the values of a,b,c in true.vals' rows
  true.logphi <- runif(ngenes, -6, -1) #values for log(g*phi)
  est.logphi <- rnorm(length(true.logphi), true.logphi, phi.sd)
  est.logphig <- rnorm(length(true.logphi), a+(b+1)*true.logphi, sqrt(c))
  recovered.vals <- optim(c(1,1,1), fn=LikelihoodForOptim, logphi.vector=est.logphi, logphig.vector=est.logphig)$par
  recovered.vals[3] <- exp(recovered.vals[3]) #optimize in log space, convert back
  est.vals <- rbind(est.vals, recovered.vals)  
  print(c(true.vals[rep,], est.vals[rep,]))
}
par(mfcol=c(1,3))
for (trait in sequence(3)) {
  true<-true.vals[,trait]
  est<-est.vals[,trait]
  plot(true, est, bty="n")
  
  fit <- lm(est ~ true -1)
  abline(fit)
  abline(a=0, b=1, lty="dotted", col="red")
  print(fit)
  legend("topleft", bty="n", legend=paste("R^2 = ", format(summary(fit)$adj.r.squared, digits=4)))
}
