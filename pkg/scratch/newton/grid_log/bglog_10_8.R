rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-7),1:106)) 
 save(res,file="bglog_10_8.RData",compress=TRUE) 
