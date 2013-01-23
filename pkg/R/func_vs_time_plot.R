## All *.RData files are stored on lab9
##################################################################################################
## read in the MLE of parameters for all genes under max rule, s_max stores the g values
##################################################################################################
maxdir <- "~/Dropbox/CollabChaiGilchristOMeara/GrantFigure2/rokas_max/"
res_max <- vector("list",length=106)
l <- 106
for(genect in 1:l){
  filename = paste(maxdir, "gene",genect,"_s_weight.RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  load(filename)
  res_max[[genect]] <- res_op
}

s_max <- sapply(1:106,function(x) res_max[[x]]$s)
rm(list=setdiff(ls(),"s_max")) #only save the object s_max

source("~/Documents/MyDocuments/Active/proteinevoutk20/pkg/R/main.R")
##################################################################################################
## Plot functionality along the simulated evolution paths, starting from sequences
## drawn from both our model and WAG model
## gene_num: order of gene; count: number of simulations to plot; margin: control ylim
## random: should the simulations be drawn from the 100 randomly or the first several?
##################################################################################################
plot.func <- function(gene_num,count=3,margin=0.1,random=TRUE){
  filename <- paste("~/Dropbox/CollabChaiGilchristOMeara/GrantFigure2/Wagprunetree/wgene",gene_num,".RData",sep="") #RData file with simulated data
  load(filename)
  obs.fty <- Ftny_protein(protein=datanum[6,],protein_op=opaa,s=s,DisMat=dismat) #functionality of observed sequence
  ## find the limits for y-axis
  fty_lim <- rbind(range(wsimWag[[1]]$fty),range(wsim[[1]]$fty),range(sim[[1]]$fty),range(simWag[[1]]$fty),rep(obs.fty,2)) ##ranges of functionalities for several typical simulations
  ylim <- c(min(fty_lim[,1]),max(fty_lim[,2])) #lower and upper bound
  ylim <- c(ylim[1]-ylim[1]*margin,ylim[2]+ylim[2]*margin) #extend them with a little margin
  
  if(random)
    index <- sample(100,count)
  else
    index <- 1:count
  ##start a plot
  #plot(sim[[index[1]]]$ftyfun,do.points=FALSE,xaxs="i",col="red",xlab="time",ylab="functionality",ylim=c(0,1),xlim=c(0,brlen),
  #     main=paste("gene",gene_num,",  g*Phi=",round(s_max[gene_num],2),sep=""), bty="n") 
  plot(sim[[index[1]]]$ftyfun,do.points=FALSE,xaxs="i",col="red",xlab="Time",ylab="Functionality",ylim=c(0,1),xlim=c(0,brlen),
       main="Functionality vs time", bty="n",xaxt="n") 
  abline(h=obs.fty,col="black",lwd=4, lty="dotted")
  for(i in index){
    plot(sim[[i]]$ftyfun,do.points=FALSE,xaxs="i",col=rgb(1,0,0,.4),add=TRUE)
    plot(simWag[[i]]$ftyfun,do.points=FALSE,xaxs="i",col=rgb(0,0,1,.4),add=TRUE)
    plot(wsim[[i]]$ftyfun,do.points=FALSE,xaxs="i",col=rgb(1,0,0,.4),add=TRUE)
    plot(wsimWag[[i]]$ftyfun,do.points=FALSE,xaxs="i",col=rgb(0,0,1,.4),add=TRUE)
  }
 # axis(1, at=seq(from=0, to=1, length.out=5)*brlen, labels=round(seq(from=0, to=1, length.out=5),2))
  axis(1, at=seq(from=0, to=1, length.out=5)*brlen, labels=round(seq(from=0, to=1, length.out=5)*brlen, 2))
  #axis(1)
}
##call of function: (gene72 has the highest g*Phi, gene83 has the lowest)
## plot.func(72,count=10,margin=0.01,random=T)
plot.func(72,count=10,margin=0.01,random=T)