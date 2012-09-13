source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_9_500_0.7.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[9]],sim[1:14,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_9_500_0.7.RData",compress=TRUE) 
