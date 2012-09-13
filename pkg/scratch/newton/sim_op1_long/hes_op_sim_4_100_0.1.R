source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_4_100_0.1.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[4]],sim[1:8,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_4_100_0.1.RData",compress=TRUE) 
