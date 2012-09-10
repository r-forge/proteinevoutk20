rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-6,-14),1:106)) 
 save(res,file="bglog_9_1.RData",compress=TRUE) 
