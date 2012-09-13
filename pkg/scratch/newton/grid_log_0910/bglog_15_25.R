rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.2,-3.16),1:106)) 
 save(res,file="bglog_15_25.RData",compress=TRUE) 
