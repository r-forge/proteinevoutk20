rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.6,-7),1:106)) 
 save(res,file="bglog_3_1.RData",compress=TRUE) 
