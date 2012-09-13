rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.199999999999999,-3.48),1:106)) 
 save(res,file="bglog_25_23.RData",compress=TRUE) 
