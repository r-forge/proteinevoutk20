rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,0),1:106)) 
 save(res,file="bglog_10_15.RData",compress=TRUE) 
