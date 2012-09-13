rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-5.24),1:106)) 
 save(res,file="bglog_5_12.RData",compress=TRUE) 
