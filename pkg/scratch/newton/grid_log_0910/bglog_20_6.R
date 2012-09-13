rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.2,-6.2),1:106)) 
 save(res,file="bglog_20_6.RData",compress=TRUE) 
