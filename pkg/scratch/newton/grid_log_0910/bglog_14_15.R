rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.4,-4.76),1:106)) 
 save(res,file="bglog_14_15.RData",compress=TRUE) 
