rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.8,-6.04),1:106)) 
 save(res,file="bglog_12_7.RData",compress=TRUE) 
