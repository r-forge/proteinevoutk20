######################################
##  need alpha, beta and gamma values to run
######################################
## load RData file that includes "data_res"
load("collect_data_log.RData")
##numer of sites in the simulation 
numsites <- c(100,300,500,700,1000)
##s values used in simulation                                                      
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
##tips of trees used in simulation                                                 
tip <- c(4,4,8,8,12,12)

## get the mle for the parameter with index paraind
get.mle <- function(treeind,siteind,sind,paraind){
  datares[[treeind,siteind,sind]]$mle[paraind]
}
vec.mle <- Vectorize(get.mle,vectorize.args=c("treeind","siteind"))
## get the eigenvalues of hessian for a simulation
get.hes.eigenval <- function(treeind,siteind,sind){
  datares[[treeind,siteind,sind]]$eigenval
}
vec.hes.eigenval <- Vectorize(get.hes.eigenval,vectorize.args=c("treeind","siteind"))
##plot the absolute error of parameter according to the number of sites
## 6 curves, one for each tree, sind specifies the index of s value,
## paraind specifies the index of parameters(1:s,2:beta,3:gamma)
plot.s <- function(sind,paraind){
  ## find mle for all the simulations, stored in an array (6*5)
  d <- outer(1:6,1:5,vec.mle,sind=sind,paraind=paraind)
  ## parameter s
  if(paraind==1){
    d <- abs(d-log(svalues[sind]))
    plot.title = paste("log(s)= ",log(svalues[sind]),sep="")
    plot.ylab = "absolute error of log(s)"
  }
  ## parameter beta/alpha
  else if(paraind==2){
    d <- abs(d-log(be))
    plot.title =  paste("s= ",svalues[sind], ", log(beta) = ", log(be), sep="")
    plot.ylab = "absolute error of log(beta)"
  }
  ## parameter gamma/alpha
  else{
    d <- abs(d-log(ga))
    plot.title =  paste("s= ",svalues[sind], ", log(gamma) = ", log(ga), sep="")
    plot.ylab = "absolute error of log(gamma)"
  }
  plot(numsites,d[1,],type="l",xlab="number of sites",ylab=plot.ylab,ylim=c(min(d),max(d)+max(d)*0.1))
  title(plot.title)
  points(numsites,d[2,],type="l",col="black")
  points(numsites,d[3,],type="l",col="blue")
  points(numsites,d[4,],type="l",col="blue")
  points(numsites,d[5,],type="l",col="red")
  points(numsites,d[6,],type="l",col="red")
  leg.txt <- c("4 tips","8 tips","12 tips")
  y.leg <- seq(from=max(d)+max(d)*0.1,by=-max(d)*0.1,length.out=3)
  col.leg <- c("black","blue","red")
  legend("topright",leg.txt,col=col.leg,text.col=col.leg,
           lty=1,bty="n",cex=0.7)
}
  
