rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-8,-7),1:106)) 
 save(res,file="bglog_7_8.RData",compress=TRUE) 
