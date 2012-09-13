rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.2,-5.08),1:106)) 
 save(res,file="bglog_10_13.RData",compress=TRUE) 
