rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.6,-3.8),1:106)) 
 save(res,file="bglog_13_21.RData",compress=TRUE) 
