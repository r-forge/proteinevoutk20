rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-10,-4),1:106)) 
 save(res,file="bglog_5_11.RData",compress=TRUE) 
