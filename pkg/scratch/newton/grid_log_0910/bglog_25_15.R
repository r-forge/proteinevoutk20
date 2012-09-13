rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.199999999999999,-4.76),1:106)) 
 save(res,file="bglog_25_15.RData",compress=TRUE) 
