rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.8,-5.4),1:106)) 
 save(res,file="bglog_17_11.RData",compress=TRUE) 
