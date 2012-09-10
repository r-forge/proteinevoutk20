rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-8),1:106)) 
 save(res,file="bglog_14_7.RData",compress=TRUE) 
