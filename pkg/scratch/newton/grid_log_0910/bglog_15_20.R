rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.2,-3.96),1:106)) 
 save(res,file="bglog_15_20.RData",compress=TRUE) 
