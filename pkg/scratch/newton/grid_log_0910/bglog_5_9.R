rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-5.72),1:106)) 
 save(res,file="bglog_5_9.RData",compress=TRUE) 
