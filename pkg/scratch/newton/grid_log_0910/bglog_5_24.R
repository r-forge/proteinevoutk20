rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-3.32),1:106)) 
 save(res,file="bglog_5_24.RData",compress=TRUE) 
