rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,0),1:106)) 
 save(res,file="bglog_12_15.RData",compress=TRUE) 
