rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.6,-6.84),1:106)) 
 save(res,file="bglog_13_2.RData",compress=TRUE) 
