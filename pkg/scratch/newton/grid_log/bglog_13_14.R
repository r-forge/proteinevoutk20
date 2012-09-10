rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2,-1),1:106)) 
 save(res,file="bglog_13_14.RData",compress=TRUE) 
