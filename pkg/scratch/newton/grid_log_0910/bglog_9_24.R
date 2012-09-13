rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.4,-3.32),1:106)) 
 save(res,file="bglog_9_24.RData",compress=TRUE) 
