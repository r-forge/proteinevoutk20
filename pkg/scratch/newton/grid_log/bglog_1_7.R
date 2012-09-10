rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-14,-8),1:106)) 
 save(res,file="bglog_1_7.RData",compress=TRUE) 
