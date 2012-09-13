rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.8,-5.72),1:106)) 
 save(res,file="bglog_17_9.RData",compress=TRUE) 
