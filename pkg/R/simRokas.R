load("~/proteinevoutk20/pkg/scratch/newton/rokas_mle_newmat/rokasMLE.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
index <- attr(res_op$data,"index")
opaa <- res_op$ll$opaa
protein_op <- opaa[index]
tree <- res_op$tree
s <- res_op$s
beta <- res_op$GMweights[2]
gamma <- res_op$GMweights[3]
bfaa <- res_op$bfaa
GTRvec <- res_op$Q
sim <- simTree(tree,protein_op,s,GTRvec,alpha=al,beta=beta,gamma=gamma,bfaa=bfaa)
save.image(file="simRokas.RData",compress=TRUE)