rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.6,-7),1:106)) 
 save(res,file="bglog_23_1.RData",compress=TRUE) 
