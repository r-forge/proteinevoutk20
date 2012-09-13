rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.4,-3.8),1:106)) 
 save(res,file="bglog_9_21.RData",compress=TRUE) 
