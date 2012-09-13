source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_14_500_1.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[14]],sim[1:18,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_14_500_1.RData",compress=TRUE) 
