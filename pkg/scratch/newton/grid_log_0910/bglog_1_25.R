rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-3.16),1:106)) 
 save(res,file="bglog_1_25.RData",compress=TRUE) 
