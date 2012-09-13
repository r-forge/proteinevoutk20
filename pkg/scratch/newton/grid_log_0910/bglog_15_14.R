rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.2,-4.92),1:106)) 
 save(res,file="bglog_15_14.RData",compress=TRUE) 
