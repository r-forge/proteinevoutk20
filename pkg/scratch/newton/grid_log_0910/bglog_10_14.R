rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.2,-4.92),1:106)) 
 save(res,file="bglog_10_14.RData",compress=TRUE) 
