load("~/proteinevoutk20/pkg/RData/sim2simple_1_100_0.1.RData")
source("~/proteinevoutk20/pkg/R/pphyproevo.R")
system.time(mle <- MLE_GTR_sw(c(0.1,0.01,0.001),rep(0,3),rep(1,3),trees[[1]],sim[1:4,],m=20,al,mumat,protein_op=op_seq,root=start_seq,trace=1))
save.image(file="sim2simpleMle_1_100_0.1.RData",compress=TRUE)
