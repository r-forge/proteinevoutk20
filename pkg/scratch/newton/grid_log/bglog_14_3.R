rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-12),1:106)) 
 save(res,file="bglog_14_3.RData",compress=TRUE) 
