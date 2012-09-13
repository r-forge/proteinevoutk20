rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.199999999999999,-5.72),1:106)) 
 save(res,file="bglog_25_9.RData",compress=TRUE) 
