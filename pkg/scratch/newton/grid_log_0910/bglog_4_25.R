rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.4,-3.16),1:106)) 
 save(res,file="bglog_4_25.RData",compress=TRUE) 
