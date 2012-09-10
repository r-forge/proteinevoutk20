rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-12,-6),1:106)) 
 save(res,file="bglog_3_9.RData",compress=TRUE) 
