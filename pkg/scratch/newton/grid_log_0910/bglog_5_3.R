rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-6.68),1:106)) 
 save(res,file="bglog_5_3.RData",compress=TRUE) 
