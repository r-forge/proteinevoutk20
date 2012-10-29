FtnyVec <- function(op,vec,s=0.1,DisMat=GM){
  sapply(vec,function(x) Ftny_protein(x,op,s,DisMat))
}

sim <- simulation(1,2,100,0.1,GM,MUMAT)
plot(sim,type="l",xlab="evolving path",ylab="Functionality",main="A simulated path when optimal aa = R")
abline(h=c(1,Ftny2[16]),col=c("red","blue"))