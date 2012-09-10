rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2,-10),1:106)) 
 save(res,file="bglog_13_5.RData",compress=TRUE) 
