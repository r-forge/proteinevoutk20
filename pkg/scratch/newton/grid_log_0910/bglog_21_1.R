rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-7),1:106)) 
 save(res,file="bglog_21_1.RData",compress=TRUE) 
