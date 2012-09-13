rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.6,-3.64),1:106)) 
 save(res,file="bglog_13_22.RData",compress=TRUE) 
