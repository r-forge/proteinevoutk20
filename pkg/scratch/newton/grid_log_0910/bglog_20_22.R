rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.2,-3.64),1:106)) 
 save(res,file="bglog_20_22.RData",compress=TRUE) 
