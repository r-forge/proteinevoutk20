rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2,-6.68),1:106)) 
 save(res,file="bglog_16_3.RData",compress=TRUE) 
