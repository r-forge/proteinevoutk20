rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2,-5.4),1:106)) 
 save(res,file="bglog_16_11.RData",compress=TRUE) 
