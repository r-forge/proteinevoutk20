num.tip <- 32
#setwd(paste("~/BackupProEvo/Newton/balBobyqa_eqm/tip",num.tip,sep=""))
load(paste("~/BackupProEvo/Newton/balBobyqa_eqm/tip",num.tip,"/",num.tip,"tip3000char.RData",sep=""))
opaa <- opaa[index]
dif <- which(opaa!=opaa1)
opaa.dif <- opaa[dif]
opaa1.dif <- opaa1[dif]

dif.data.AA <- matrix(AA[datanum[,dif]],nrow=num.tip)
rownames(dif.data.AA) <- rownames(datanum)
dif.data <- phyDat(dif.data.AA,type="AA")
Qall <- res.op$Qall
tree <- res.op$tree
dif.ind <- attr(dif.data,"index")
opaa_cmp <- function(y)
  sapply(c(opaa[dif[y]],opaa1[dif[y]]),
         function(x) exp(as.vector(ll3m(dif.data,tree,Q=Qall[[x]])$re))[dif.ind[y]])
#################################################################
opaas <- cbind(opaa1[dif],opaa[dif]) #only the cases where inferred opaa is different from true values
#opaas.unique <- cbind(opaa1,opaa)
opaas.unique <- unique.array(opaas) #unique rows
opaas.unique <- opaas.unique[order(opaas.unique[,1]),] #order by the true values
counts <- NULL #counts of each mismatch pair
for(i in 1:dim(opaas.unique)[1]){
  count <- length(intersect(which(opaas[,1]==opaas.unique[i,1]),which(opaas[,2]==opaas.unique[i,2])))
  counts <- c(counts,count)
}
#par(mfrow=c(2,1),oma = c(0, 0, 2, 0))
#bubble plot 
symbols(x=opaas.unique[,1],y=opaas.unique[,2],circles=counts,inches=1/4,bg="red",
        xlim=c(1,20),ylim=c(1,20),main=paste(num.tip,"taxa"),
        xlab="true opaa",ylab="inferred opaa",xaxt='n',yaxt='n')
# axis(1,at=1:20,labels=AA,tick=FALSE,lwd=0.5,cex.axis=0.6)
# axis(2,at=1:20,labels=AA,tick=FALSE,lwd=0.5,cex.axis=0.6)
# abline(v=1:20,h=1:20)
#################################################################
opaas.all <- cbind(opaa1,opaa)
opaas.all.unique <- unique.array(opaas.all) #unique rows
opaas.all.unique <- opaas.all.unique[order(opaas.all.unique[,1]),] #order by the true values
all.counts <- NULL #counts of each mismatch pair
for(i in 1:dim(opaas.all.unique)[1]){
  count <- length(intersect(which(opaas.all[,1]==opaas.all.unique[i,1]),which(opaas.all[,2]==opaas.all.unique[i,2])))
  all.counts <- c(all.counts,count)
}
#bubble plot 
# symbols(x=opaas.all.unique[,1],y=opaas.all.unique[,2],circles=all.counts,inches=1/4,bg="blue",
#         main="Plot of opaas",xlab="start opaa",ylab="inferred opaa",
#         xlim=c(1,20),ylim=c(1,20),xaxt='n',yaxt='n',add=TRUE)
symbols(x=opaas.all.unique[,1],y=opaas.all.unique[,2],circles=all.counts,inches=1/4,bg="blue",
        xlim=c(1,20),ylim=c(1,20),xaxt='n',yaxt='n',add=TRUE)
axis(1,at=1:20,labels=AA,tick=FALSE,lwd=0.5,cex.axis=0.6)
axis(2,at=1:20,labels=AA,tick=FALSE,lwd=0.5,cex.axis=0.6)
abline(v=1:20,h=1:20)
#mtext(paste(num.tip,"taxa"),outer=TRUE,cex=1.2,lwd=1.3)


###########################################################################
ftny_vec <- Vectorize(FUN="Ftny_protein",vectorize.args=c("protein","protein_op"),SIMPLIFY=TRUE)
###########################################################################

## plot the simulation path of one amino acid
sim8 <- simulation(protein=8,protein_op=8,t=100,s=s,DisMat=dismat,MuMat=mumat)
plot(sim8$path[,1]~sim8$path[,2],type="b",pch=20,xlab="time",ylab="amino acid",yaxt='n',ylim=c(1,20),
     main=paste("opaa =", AA[8]))
axis(2,at=1:20,labels=AA,lwd=0.5,cex.axis=0.6)
abline(h=unique(sim8$path[,1]),lty=2)