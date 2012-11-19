load("~/proteinevoutk20/pkg/scratch/newton/rokas_max/gene1_s_weight.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
#estimates of parameters from gene1 data
beta = res_op$GMweights[2]
gamma = res_op$GMweights[3]
s = res_op$s
tree = res_op$tree
index <- attr(res_op$data,"index")
opaa <- res_op$ll$opaa
protein_op <- opaa[index]
GTRvec = res_op$Q

simgene1 <- simTree(tree,protein_op,s,GTRvec,alpha=al,beta=beta,gamma=gamma,
                    rootseq=sample(20,length(protein_op),replace=T))
simgene1phy = phyDat(simgene1$data,type="AA")
str(simgene1phy)
rootseq = simgene1$rootseq
rootdis = sapply(1:20,function(x) rootseq[opaa==x])
