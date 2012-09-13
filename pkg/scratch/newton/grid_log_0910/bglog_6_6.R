rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4,-6.2),1:106)) 
 save(res,file="bglog_6_6.RData",compress=TRUE) 
