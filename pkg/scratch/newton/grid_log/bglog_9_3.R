rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-6,-12),1:106)) 
 save(res,file="bglog_9_3.RData",compress=TRUE) 
