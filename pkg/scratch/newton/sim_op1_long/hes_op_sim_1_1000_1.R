source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_1_1000_1.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[1]],sim[1:4,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_1_1000_1.RData",compress=TRUE) 
