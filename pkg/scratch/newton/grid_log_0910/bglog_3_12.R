rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.6,-5.24),1:106)) 
 save(res,file="bglog_3_12.RData",compress=TRUE) 
