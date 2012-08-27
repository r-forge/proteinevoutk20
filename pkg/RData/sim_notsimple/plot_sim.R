## load RData file that includes "data_res"
load("collect_data_log.RData")
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
  d <- outer(1:6,1:5,vec.mle,sind=sind,paraind=paraind)
  if(paraind==1){
    d <- abs(d-log(svalues[sind]))
    plot.title = paste("s= ",svalues[sind],sep="")
    plot.ylab = "absolute error of s"
  }
  else if(paraind==2){
    d <- abs(d-log(be))
    plot.title =  paste("s= ",svalues[sind], ", log(beta) = ", log(be), sep="")
    plot.ylab = "absolute error of beta"
  }
  else{
    d <- abs(d-log(ga))
    plot.title =  paste("s= ",svalues[sind], ", log(gamma) = ", log(ga), sep="")
    plot.ylab = "absolute error of gamma"
  }
  plot(numsites,d[1,],type="l",xlab="number of sites",ylab=plot.ylab,ylim=c(0,max(d)+max(d)*0.1))
  title(plot.title)
  points(numsites,d[2,],type="l",col="black")
  points(numsites,d[3,],type="l",col="blue")
  points(numsites,d[4,],type="l",col="blue")
  points(numsites,d[5,],type="l",col="red")
  points(numsites,d[6,],type="l",col="red")
  leg.txt <- c("4 tips","8 tips","12 tips")
  y.leg <- seq(from=max(d)+max(d)*0.1,by=-max(d)*0.1,length.out=3)
  col.leg <- c("black","blue","red")
  for(i in 1:3)
    legend(100,y.leg[i],leg.txt[i],col=col.leg[i],text.col=col.leg[i],
           lty=1,bty="n",cex=0.7)
}
  
