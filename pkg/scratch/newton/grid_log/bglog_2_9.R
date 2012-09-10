rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-13,-6),1:106)) 
 save(res,file="bglog_2_9.RData",compress=TRUE) 
