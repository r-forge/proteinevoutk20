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

sim1 <- simTree(tree,rep(1,400),s,GTRvec,al,be,ga,rootseq=rep(1:20,20))
sim1phy <- phyDat(sim1$data,type="AA")
str(sim1phy)

#find the most frequent patterns from data of phyDat format 
MostFreq <- function(data,count=10){
  weight = attr(data,"weight")
  datamat <- matrix(unlist(data),nrow=length(data),byrow=T)
  Index = head(order(weight,decreasing=T),n=count)
  return(list(mat=datamat[,Index],wt=weight[Index]))
}