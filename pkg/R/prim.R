source("~/proteinevoutk20/pkg/R/main.R")
prim <- read.nexus.data("~/proteinevoutk20/pkg/Data/primates.nex")
for(i in 1:length(prim)){
  prim[[i]] <- prim[[i]][c(2:457,660:896)] #coding region
  }
prim.AA <- lapply(prim,translate,numcode=2) #mtDNA
prim.phy <- phyDat(prim.AA,type="AA") #in phyDat format
tre <- read.nexus("~/proteinevoutk20/pkg/Data/primates.tre") #read the tree
rtre <- (root(tre,1:2,resolve.root=T)) # root the phylogeny
prim.vec <- c(3.56175, 43.70565, 3.02842, 2.34962, 34.63963, 1.00000)
prim.res <- mllm1(prim.phy,rtre,s=1,beta=be,gamma=ga,Q=rep(1,6))
optim.prim <- optim.mllm1(prim.res,optQ=T,optBranch=T,optsWeight=T,control = list(epsilon=1e-08,hmaxit=20,htrace=TRUE,print_level=1,maxeval="200"))
save(prim.phy,optim.prim,file="prim.RData")

         
