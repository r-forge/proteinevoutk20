rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.6,-5.4),1:106)) 
 save(res,file="bglog_13_11.RData",compress=TRUE) 
