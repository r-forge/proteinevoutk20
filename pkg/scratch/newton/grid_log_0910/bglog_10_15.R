rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.2,-4.76),1:106)) 
 save(res,file="bglog_10_15.RData",compress=TRUE) 
