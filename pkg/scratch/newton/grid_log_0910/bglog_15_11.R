rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.2,-5.4),1:106)) 
 save(res,file="bglog_15_11.RData",compress=TRUE) 
