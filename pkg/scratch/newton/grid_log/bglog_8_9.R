rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-7,-6),1:106)) 
 save(res,file="bglog_8_9.RData",compress=TRUE) 
