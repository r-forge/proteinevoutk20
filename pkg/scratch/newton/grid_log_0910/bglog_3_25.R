rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.6,-3.16),1:106)) 
 save(res,file="bglog_3_25.RData",compress=TRUE) 
