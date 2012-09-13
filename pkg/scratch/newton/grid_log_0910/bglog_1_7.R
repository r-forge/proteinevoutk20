rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-6.04),1:106)) 
 save(res,file="bglog_1_7.RData",compress=TRUE) 
