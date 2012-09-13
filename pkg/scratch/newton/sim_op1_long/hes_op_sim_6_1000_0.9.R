source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_6_1000_0.9.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[6]],sim[1:10,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_6_1000_0.9.RData",compress=TRUE) 
