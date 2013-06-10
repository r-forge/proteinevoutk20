load("~/proteinevoutk20/pkg/Data/Rokas/rokas.RData")
## including AA matrix called "rokas" and range for each gene called "AArange"
############################################################################
# read in Phi values for 106 (and more) genes in budding yeast
# rokasPhi <- read.csv("~/proteinevoutk20/pkg/Data/Rokas/genePhi.csv",header=TRUE)
# phi <- read.csv("~/proteinevoutk20/pkg/Data/phi.csv")
# # convert class of data frame from factor to character
# phi <- data.frame(lapply(phi,as.character),stringsAsFactors=FALSE)
# rokasPhi <- data.frame(lapply(rokasPhi,as.character),stringsAsFactors=FALSE)
# names(rokasPhi)[2] <- "Ingolia" #Ingolia Phi
# rokasPhi$Beyer <- NA #Beyer Phi
# rokasPhi$SEMPPR <- NA #SEMPPR Phi
# rokasPhi$Ing <- NA #Ingolia Phi from data that including all 3 columns
# for(i in 1:106){
#   if(rokasPhi$gene[i] %in% phi$Orf){ #If this gene is inclued in the 3 column data
#     rokasPhi$Beyer[i] <- phi$Beyer_Phi[which(phi$Orf == rokasPhi$gene[i])]
#     rokasPhi$SEMPPR[i] <- phi$SEMPPR_Mean_Phi[which(phi$Orf == rokasPhi$gene[i])]
#     rokasPhi$Ing[i] <- phi$Ingolia_Phi[which(phi$Orf == rokasPhi$gene[i])]
#   }
# }
# #rm(phi)
# rokasPhi$Ingolia <- as.numeric(rokasPhi$Ingolia) #convert character class to numeric class
# rokasPhi$Beyer <- as.numeric(rokasPhi$Beyer)
# rokasPhi$SEMPPR <- as.numeric(rokasPhi$SEMPPR)
# rokasPhi$Ing <- as.numeric(rokasPhi$Ing)

############################################################################
## analysis of mtDNA data of primates on 12 species, length = 898 including noncoding region
# prim <- read.nexus.data("~/Garli-2.0-anyOSX/example/basic/primates.nex")
# for(i in 1:length(prim)){
#   prim[[i]] <- prim[[i]][c(2:457,660:896)] #coding region
# }
# prim.AA <- lapply(prim,translate,numcode=2) #mtDNA
# prim.phy <- phyDat(prim.AA,type="AA") #in phyDat format
# tre <- read.nexus("~/proteinevoutk20/pkg/Data/primates.tre") #read the tree
# rtre <- (root(tre,1:2,resolve.root=T)) # root the phylogeny
# prim.vec <- c(3.56175, 43.70565, 3.02842, 2.34962, 34.63963, 1.00000)
# prim.res <- mllm1(prim.phy,rtre,s=1,beta=be,gamma=ga,Q=rep(1,6))
# optim.prim <- optim.mllm1(prim.res,optQ=T,optBranch=T,optsWeight=T,control = list(epsilon=1e-08,hmaxit=10,htrace=TRUE,print_level=1,maxeval="300"))

