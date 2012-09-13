rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.8,-7),1:106)) 
 save(res,file="bglog_17_1.RData",compress=TRUE) 
