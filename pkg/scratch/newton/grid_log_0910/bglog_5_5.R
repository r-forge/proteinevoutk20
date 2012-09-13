rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-6.36),1:106)) 
 save(res,file="bglog_5_5.RData",compress=TRUE) 
