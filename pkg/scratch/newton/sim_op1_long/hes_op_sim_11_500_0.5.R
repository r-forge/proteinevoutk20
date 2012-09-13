source("~/proteinevoutk20/pkg/R/hessian.R") 
 load("op_simu_11_500_0.5.RData") 
 hes_op <- find_hessian_log(mle_log_op$par,trees[[11]],sim[1:16,],al,mumat,20,protein_op=op_seq,root=start_seq) 
 save.image(file="hes_op_simu_11_500_0.5.RData",compress=TRUE) 
