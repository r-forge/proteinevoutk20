##   Read in 106 gene data  from Rokas's ##
#source("~/proteinevoutk20/pkg/R/main.R")
l <- 106
ROKAS_DATA <- vector("list",length=l)
##data stores the data of all genes in rokas's data
for(i in 1:l){
  ROKAS_DATA[[i]] <- conv(paste(datadir,"gene",i,".fasta",sep=""),type="phyDat")
}
rm(l)
rokasPhi <- read.csv("~/proteinevoutk20/pkg/Data/genePhi.csv",header=TRUE)
rokasdata = read.phyDat("~/proteinevoutk20/pkg/Data/rokasAA",type="AA")

phi <- read.csv("Desktop/phi.csv")
# convert class of data frame from factor to character
phi <- data.frame(lapply(phi,as.character),stringsAsFactors=FALSE)
rokasPhi <- data.frame(lapply(rokasPhi,as.character),stringsAsFactors=FALSE)
names(rokasPhi)[2] <- "Ingolia" #Ingolia Phi
rokasPhi$Beyer <- NA #Beyer Phi
rokasPhi$SEMPPR <- NA #SEMPPR Phi
for(i in 1:106){
  if(rokasPhi$gene[i] %in% phi$Orf){ #If this gene is inclued in the 3 column data
    rokasPhi$Beyer[i] <- phi$Beyer_Phi[which(phi$Orf == rokasPhi$gene[i])]
    rokasPhi$SEMPPR[i] <- phi$SEMPPR_Mean_Phi[which(phi$Orf == rokasPhi$gene[i])]
  }
}
rm(phi)
rokasPhi$Ingolia <- as.numeric(rokasPhi$Ingolia) #convert character class to numeric class
rokasPhi$Beyer <- as.numeric(rokasPhi$Beyer)
rokasPhi$SEMPPR <- as.numeric(rokasPhi$SEMPPR)