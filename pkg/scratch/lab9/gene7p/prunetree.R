load("gene7p.RData")
source("~/proteinevoutk20/pkg/R/main.R")
## the objects for the pruned data/tree has suffix "p" indicating that it's for the pruned
## use the paraemters found with 7 tips, find the optimal branch lengths with 8 tips
s = res_op$s
beta= res_op$GMweights[2]
gamma = res_op$GMweights[3]
Q = res_op$Q
treep = res_op$tree
opaa = res_op$ll$opaa #with data that only have 7 tips, cannot be used in the mle search for branch 
                      #length below  
gene7 <- conv("~/proteinevoutk20/pkg/Data/gene7.fasta",type="phyDat")
gene7p <- conv("~/proteinevoutk20/pkg/scratch/lab9/gene7p/gene7p.fasta",type="phyDat")

tree <- ROKAS_TREE
br <- optim.br(gene7,tree,maxeval="3000",print_level=1,s=s,beta=beta,gamma=gamma,Q=Q)
br
save.image("gene7_8tips.RData)
