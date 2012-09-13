rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.6,-4.6),1:106)) 
 save(res,file="bglog_8_16.RData",compress=TRUE) 
