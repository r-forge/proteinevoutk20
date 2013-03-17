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
#rokasPhi$Ing <- rokasPhi$Ingolia/10000
############################################################################
phi <- rokasPhi
phi$Ing <- phi$Ingolia/10000 # less NA's in Ing column
phi$Ingolia <- NULL #drop the first column, which is proportional to the last column
phi$gene <- NULL #drop gene names
na.col <- sort(unique(arrayInd(which(is.na(phi)),.dim=c(106,3))[,1])) #columns with NA in it
phi <- phi[-na.col,]
phi.mat <- data.matrix(phi) #convert to matrix
gphi <- gphi[-na.col] # phi and gphi without NA's 
lgphi <- log(gphi)

phidata <- log(phi)
phidata$lgphi <- lgphi #put all log(phi) values and log(gphi) in the same data frame
allphi <- c(data.matrix(phidata[,1:3]))
############################################################################
## same data frame with log(phi) values centered at 0, this does not affect the 
## coefficients in the linear regression, but does affect the intercept
## new data frame is called cphidata, also with a column for log(gphi) values
############################################################################
phi.log <- log(phi) # log phi
phi.log.center <- phi.log
for(i in 1:3){
  phi.log.center[,i] <- phi.log[,i] - mean(phi.log[,i],na.rm=TRUE) #centered log phi values
}
#phi.log.mean <- colMeans(phi.log.center,na.rm=TRUE)
clogphi <- phi.log.center #rename the centered log phi values as cphi -- shorter name
cphidata <- clogphi
cphidata$lgphi <- lgphi
allcphi <- c(data.matrix(cphidata[,1:3]))
############################################################################
## analysis procedure on data frame structured like the cphidata or phidata
## the data frame should have the same structure as phidata,
## 1st col: Beyer; 2nd col: SEMPPR; 3rd col: Ing; 4th col: log(gphi)
############################################################################
phidata.lm <- function(phidata){
  lm.phi <- lm(lgphi~Beyer+SEMPPR+Ing,data=phidata) #linear regression on 3 measurements
  #print(summary(lm.phi)) #print regression info
  predict.gphi <- lm.phi$fitted.values # est(gphi)
  weighted.phi <- data.matrix(phidata[,1:3]) %*% lm.phi$coefficients[-1] #weighted phi
  est.gval <- phidata$lgphi - weighted.phi
  #est.gval <- lm.phi$coefficients[1]+weighted.phi #est(g)
  gval <- predict.gphi-weighted.phi
  #gval <- lm.phi$residuals + lm.phi$coefficients[1] #g values
  opar <- par(no.readonly=TRUE) #original graphical parameter
  par(mfrow=c(2,2)) # new graphical parameter
  ##Plot of g vs. phi values
  plot(est.gval~weighted.phi,main="g vs weighted phi/log scale",xlab="weighted phi",ylab="g")
  plot(est.gval~phidata$Beyer,main="g vs Beyer phi/log scale",xlab="Beyer phi",ylab="g")
  plot(est.gval~phidata$SEMPPR,main="g vs SEMPPR phi/log scale",xlab="SEMPPR phi",ylab="g")
  plot(est.gval~phidata$Ing,main="g vs Ingolia phi/log scale",xlab="Ingolia phi",ylab="g")
  par(opar) #reset to original
  return(lm.phi)
}

lm.plot <- function(lm,phi){
  predict.gphi <- lm$fitted.values # est(gphi)
  #weighted.phi <- data.matrix(phidata[,1:3]) %*% lm.phi$coefficients[-1] #weighted phi
  est.gval <- lm$coefficients[1]+(lm$coefficients[2]-1)*lm$model[[2]] #est(g): a + b*log(phi)
  plot((predict.gphi-phi)~phi,main="est(gphi)/obs(phi) vs. obs(phi) on log scale",xlab="obs(phi)",ylab="est(gphi)/obs(phi)")
  plot(est.gval~phi,main="est(g) vs. obs(phi) on log scale",xlab="obs(phi)",ylab="est(g)")
}