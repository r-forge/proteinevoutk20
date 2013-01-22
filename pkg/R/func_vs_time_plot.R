## All *.RData files are stored on lab9
##################################################################################################
## read in the MLE of parameters for all genes under max rule, s_max stores the g values
##################################################################################################
maxdir <- "~/BackupProEvo/Newton/rokas_max/"
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

source("~/proteinevoutk20/pkg/R/main.R")
##################################################################################################
## Plot functionality along the simulated evolution paths, starting from sequences
## drawn from both our model and WAG model
## gene_num: order of gene; count: number of simulations to plot; margin: control ylim
## random: should the simulations be drawn from the 100 randomly or the first several?
##################################################################################################
plot.func <- function(gene_num,count=3,margin=0.1,random=TRUE){
  filename <- paste("~/BackupProEvo/Lab9/Wagprunetree/wgene",gene_num,".RData",sep="") #RData file with simulated data
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
  plot(sim[[index[1]]]$ftyfun,do.points=FALSE,xaxs="i",col="red",xlab="time",ylab="functionality",ylim=ylim,xlim=c(0,brlen),
       main=paste("gene",gene_num,",  g*Phi=",round(s_max[gene_num],2),sep="")) 
  abline(h=obs.fty,col="green",lwd=2) 
  for(i in index){
    plot(sim[[i]]$ftyfun,do.points=FALSE,xaxs="i",col="red",add=TRUE)
    plot(simWag[[i]]$ftyfun,do.points=FALSE,xaxs="i",col="blue",add=TRUE)
    plot(wsim[[i]]$ftyfun,do.points=FALSE,xaxs="i",col="red",add=TRUE)
    plot(wsimWag[[i]]$ftyfun,do.points=FALSE,xaxs="i",col="blue",add=TRUE)
  }
}
##call of function: (gene72 has the highest g*Phi, gene83 has the lowest)
## plot.func(1,count=10,margin=0.01,random=T)