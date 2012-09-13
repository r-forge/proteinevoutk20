rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.8,-6.52),1:106)) 
 save(res,file="bglog_17_4.RData",compress=TRUE) 
