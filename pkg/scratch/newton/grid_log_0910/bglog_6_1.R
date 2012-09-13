rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4,-7),1:106)) 
 save(res,file="bglog_6_1.RData",compress=TRUE) 
