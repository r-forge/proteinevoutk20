rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.6,-3.48),1:106)) 
 save(res,file="bglog_8_23.RData",compress=TRUE) 
