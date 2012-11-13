FtnyVec <- function(op,vec,s=0.1,DisMat=GM){
  sapply(vec,function(x) Ftny_protein(x,op,s,DisMat))
}

sim <- simulation(1,2,100,0.1,GM,MUMAT)
plot(sim,type="l",xlab="evolving path",ylab="Functionality",main="A simulated path when optimal aa = R")
abline(h=c(1,Ftny2[16]),col=c("red","blue"))

modeAA <- function(data){
  return(sapply(1:dim(data)[2],function(x) Mode(data[,x])))
}
#number of distinct AA's observed at different species for 1 site
disAA <- function(data){
  return(sapply(1:dim(data)[2],function(x) length(unique(data[,x]))))
}

modeaa <- sapply(1:1000,function(x) Mode(sim1$data[,x]))
disaa <- sapply(1:1000,function(x) length(unique(sim1$data[,x])))

for(i in 1:106){
  Phivals[Phivals[,1]==genenames[i]]
}