rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.6,-3.32),1:106)) 
 save(res,file="bglog_13_24.RData",compress=TRUE) 
