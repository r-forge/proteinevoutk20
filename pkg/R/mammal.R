source("~/proteinevoutk20/pkg/R/main.R")
#mamdata <- read.nexus.data("proteinevoutk20/pkg/Data/mammals/M15taxa.nex")
#write.fasta(mam1,names(mam1),"~/proteinevoutk20/pkg/Data/mammals/mam1.fasta",nbchar=46152)
mamAA <- conv("~/proteinevoutk20/pkg/Data/mammals/mam1.fasta",range=NULL,type="phyDat")
mamtree <- read.nexus("proteinevoutk20/pkg/Data/mammals/T15taxa.nex")
mam.res <- mllm1(mamAA,mamtree,s=1,beta=be,gamma=ga,Q=rep(1,6))
mam.optim <- optim.mllm1(mam.res,optQ=T,optBranch=T,optsWeight=T,
                        control=list(epsilon=1e-08,hmaxit=30,htrace=TRUE,print_level=0,maxit=300))
save(be.optim,file="mam15taxa.RData")