##   Read in 106 gene data  from Rokas's ##
source("~/proteinevoutk20/pkg/R/main.R")
l <- 106
ROKAS_DATA <- vector("list",length=l)
##data stores the data of all genes in rokas's data
for(i in 1:l){
  ROKAS_DATA[[i]] <- conv(paste(datadir,"gene",i,".fasta",sep=""),type="phyDat")
}
rm(l)
rokasPhi <- read.csv("~/proteinevoutk20/pkg/Data/genePhi.csv",header=TRUE)
rokasdata = read.phyDat("proteinevoutk20/pkg/Data/rokasAA",type="AA")
