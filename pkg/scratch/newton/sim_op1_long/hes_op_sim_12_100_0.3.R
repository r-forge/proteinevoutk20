source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_12_100_0.3.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[12]],sim[1:16,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_12_100_0.3.RData",compress=TRUE) 
