rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.8,-4.6),1:106)) 
 save(res,file="bglog_17_16.RData",compress=TRUE) 
