######################################
##  need alpha, beta and gamma values to run
######################################
## load RData file that includes "data_res"
##load("collect_data_log.RData")
##numer of sites in the simulation 
numsites <- c(100,300,500,700,1000)
##s values used in simulation                                                      
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
##tips of trees used in simulation                                                 
tip <- c(4,4,8,8,10,10,12,12,14,14,14,14,14,14,16,16,18,18)

## get the mle for the parameter with index paraind
get.mle <- function(treeind,siteind,sind,paraind,op=TRUE){
  if(op)
    datares[[treeind,siteind,sind]]$mle_op[paraind]
  else
    datares[[treeind,siteind,sind]]$mle_nop[paraind]
}
vec.mle <- Vectorize(get.mle,vectorize.args=c("treeind","siteind"))

##plot the absolute error of parameter according to the number of sites
## 6 curves, one for each tree, sind specifies the index of s value,
## paraind specifies the index of parameters(1:s,2:beta,3:gamma)
plot.s <- function(sind,paraind,op=T){
  ## find mle for all the simulations, stored in an array (6*5)
  d <- outer(1:14,1:5,vec.mle,sind=sind,paraind=paraind,op=op)
  ## parameter s                                                                                                     
  if(paraind==1){
    d <- abs(d-log(svalues[sind]))
    plot.title = paste("log(s)= ",log(svalues[sind])," Op=", op, sep="")
    plot.ylab = "absolute error of log(s)"
  }
  ## parameter beta/alpha                                                                                            
  else if(paraind==2){
    d <- abs(d-log(be))
    plot.title =  paste("s= ",svalues[sind], ", log(beta) = ", log(be), " Op=", op,sep="")
    plot.ylab = "absolute error of log(beta)"
  }
  ## parameter gamma/alpha                                                                                           
  else{
    d <- abs(d-log(ga))
    plot.title =  paste("s= ",svalues[sind], ", log(gamma) = ", log(ga), " Op=", op,sep="")
    plot.ylab = "absolute error of log(gamma)"
  }
  
  plot(numsites,d[1,],type="l",xlab="number of sites",ylab=plot.ylab,ylim=c(min(d),min(max(d)+max(d)*0.1,2.5)))
  title(plot.title)
  points(numsites,d[2,],type="l",col="black")
  points(numsites,d[3,],type="l",col="blue")
  points(numsites,d[4,],type="l",col="blue")
  points(numsites,d[5,],type="l",col="red")
  points(numsites,d[6,],type="l",col="red")
  points(numsites,d[7,],type="l",col="yellow")
  points(numsites,d[8,],type="l",col="yellow")
  points(numsites,d[9,],type="l",col="tan")
  points(numsites,d[10,],type="l",col="tan")
  points(numsites,d[11,],type="l",col="tan")
  points(numsites,d[12,],type="l",col="tan")
  points(numsites,d[13,],type="l",col="tan")
  points(numsites,d[14,],type="l",col="tan")
#   points(numsites,d[15,],type="l",col="darkviolet")
#   points(numsites,d[16,],type="l",col="darkviolet")
#   points(numsites,d[17,],type="l",col="seagreen")
#   points(numsites,d[18,],type="l",col="seagreen")
  #leg.txt <- c("4 tips","8 tips","10 tips","12 tips","14 tips","16 tips","18tips")
  leg.txt <- c("4 tips","8 tips","10 tips","12 tips","14 tips")
  #y.leg <- seq(from=max(d)+max(d)*0.1,by=-max(d)*0.1,length.out=7)
  #col.leg <- c("black","blue","red","yellow","tan","darkviolet","seagreen")
  col.leg <- c("black","blue","red","yellow","tan")
  legend("topright",leg.txt,col=col.leg,text.col=col.leg,
           lty=1,bty="n",y.intersp=0.5,cex=0.5)
}
  
plot.s.mean <- function(sind,paraind,op=T){
  ## find mle for all the simulations, stored in an array (6*5)
  d <- outer(1:14,1:5,vec.mle,sind=sind,paraind=paraind,op=op)
  ## parameter s                                                                                                     
  if(paraind==1){
    d <- abs(d-log(svalues[sind]))
    plot.title = paste("log(s)= ",log(svalues[sind])," Op=", op, sep="")
    plot.ylab = "absolute error of log(s)"
  }
  ## parameter beta/alpha                                                                                            
  else if(paraind==2){
    d <- abs(d-log(be))
    plot.title =  paste("s= ",svalues[sind], ", log(beta) = ", log(be), " Op=", op,sep="")
    plot.ylab = "absolute error of log(beta)"
  }
  ## parameter gamma/alpha                                                                                           
  else{
    d <- abs(d-log(ga))
    plot.title =  paste("s= ",svalues[sind], ", log(gamma) = ", log(ga), " Op=", op,sep="")
    plot.ylab = "absolute error of log(gamma)"
  }
  
  plot(numsites,d[1,],type="l",xlab="number of sites",ylab=plot.ylab,ylim=c(min(d),max(d)+max(d)*0.1))
  title(plot.title)
  points(numsites,apply(d[3:4,],2,mean),type="l",col="blue")
  points(numsites,apply(d[5:6,],2,mean),type="l",col="red")
  points(numsites,apply(d[7:8,],2,mean),type="l",col="yellow")
  points(numsites,apply(d[9:10,],2,mean),type="l",col="tan")
  points(numsites,apply(d[11:12,],2,mean),type="l",col="darkviolet")
  points(numsites,apply(d[13:14,],2,mean),type="l",col="seagreen")
  leg.txt <- c("4 tips","8 tips","10 tips","12 tips","14 tips","16 tips","18tips")
  #leg.txt <- c("4 tips","8 tips","10 tips","12 tips","14 tips")

  col.leg <- c("black","blue","red","yellow","tan","darkviolet","seagreen")
  #col.leg <- c("black","blue","red","yellow","tan")
  legend("topright",leg.txt,col=col.leg,text.col=col.leg,
         lty=1,bty="n",y.intersp=0.5,cex=0.5)
}

plotall <- function(sind,avg=FALSE){
  filename = paste("s",svalues[sind],".pdf",sep="")
  pdf(file=filename)
  par(mfrow=c(3,2))
  for(k in 1:3){
    if(avg){
      plot.s.mean(sind,k)
      plot.s.mean(sind,k,F)
    }
    else{
      plot.s(sind,k)
      plot.s(sind,k,F)
    }
  }
  dev.off()
}