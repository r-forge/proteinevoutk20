rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.6,-6.04),1:106)) 
 save(res,file="bglog_3_7.RData",compress=TRUE) 
