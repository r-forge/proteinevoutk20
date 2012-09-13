rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.4,-5.4),1:106)) 
 save(res,file="bglog_9_11.RData",compress=TRUE) 
