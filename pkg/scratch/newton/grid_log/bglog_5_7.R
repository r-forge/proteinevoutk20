rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-10,-8),1:106)) 
 save(res,file="bglog_5_7.RData",compress=TRUE) 
