rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.199999999999999,-3.96),1:106)) 
 save(res,file="bglog_25_20.RData",compress=TRUE) 
