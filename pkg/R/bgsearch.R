load("~/proteinevoutk20/pkg/scratch/newton/rokas_mle_newmat/rokasMLE.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
source("~/proteinevoutk20/pkg/R/readRokas.R")
index <- attr(res_op$data,"index")
opaa <- res_op$ll$opaa
protein_op <- opaa[index]
length(protein_op)
tree <- res_op$tree
tree
s <- res_op$s
s
beta <- res_op$GMweights[2]
beta
gamma <- res_op$GMweights[3]
gamma
bfaa <- res_op$bfaa
bfaa
GTRvec <- res_op$Q
GTRvec
bg <- optim.w(beta,gamma,1:10,tree,trace=1,multicore=TRUE,Q=GTRvec)
save.image(file="bgsearch.RData",compress=TRUE)
