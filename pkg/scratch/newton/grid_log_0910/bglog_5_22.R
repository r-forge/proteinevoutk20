rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-3.64),1:106)) 
 save(res,file="bglog_5_22.RData",compress=TRUE) 
