## load RData file for gphi values
load("~/proteinevoutk20/pkg/RData/gphi.RData")
gphi <- s_max
############################################################################
## read in Phi values for 106 (and more) genes in budding yeast
rokasPhi <- read.csv("~/proteinevoutk20/pkg/Data/genePhi.csv",header=TRUE)
phi <- read.csv("~/proteinevoutk20/pkg/Data/phi.csv")
# convert class of data frame from factor to character
phi <- data.frame(lapply(phi,as.character),stringsAsFactors=FALSE)
rokasPhi <- data.frame(lapply(rokasPhi,as.character),stringsAsFactors=FALSE)
names(rokasPhi)[2] <- "Ingolia" #Ingolia Phi
rokasPhi$Beyer <- NA #Beyer Phi
rokasPhi$SEMPPR <- NA #SEMPPR Phi
rokasPhi$Ing <- NA #Ingolia Phi from data that including all 3 columns
for(i in 1:106){
  if(rokasPhi$gene[i] %in% phi$Orf){ #If this gene is inclued in the 3 column data
    rokasPhi$Beyer[i] <- phi$Beyer_Phi[which(phi$Orf == rokasPhi$gene[i])]
    rokasPhi$SEMPPR[i] <- phi$SEMPPR_Mean_Phi[which(phi$Orf == rokasPhi$gene[i])]
    rokasPhi$Ing[i] <- phi$Ingolia_Phi[which(phi$Orf == rokasPhi$gene[i])]
  }
}
#rm(phi)
rokasPhi$Ingolia <- as.numeric(rokasPhi$Ingolia) #convert character class to numeric class
rokasPhi$Beyer <- as.numeric(rokasPhi$Beyer)
rokasPhi$SEMPPR <- as.numeric(rokasPhi$SEMPPR)
rokasPhi$Ing <- as.numeric(rokasPhi$Ing)
############################################################################
phi <- rokasPhi
phi$Ingolia <- NULL #drop the first column, which is proportional to the last column
phi$gene <- NULL #drop gene names
phi.log <- log(phi) # log phi
phi.log.center <- phi.log
for(i in 1:3){
  phi.log.center[,i] <- phi.log[,i] - mean(phi.log[,i],na.rm=TRUE) #centered log phi values
}
phi.log.mean <- colMeans(phi.log.center,na.rm=TRUE)

nas <- which(is.na(phi.log.center$SEMPPR)) #genes with NA in SEMPPR phi values
gphi <- gphi[-nas] #drop those genes without phi values
phi.log.center <- phi.log.center[-nas,] #drop the same genes
cphi <- phi.log.center #rename the centered log phi values as cphi -- shorter name