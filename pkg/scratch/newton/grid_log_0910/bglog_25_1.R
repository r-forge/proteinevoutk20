rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.199999999999999,-7),1:106)) 
 save(res,file="bglog_25_1.RData",compress=TRUE) 
