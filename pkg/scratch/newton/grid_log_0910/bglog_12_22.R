rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.8,-3.64),1:106)) 
 save(res,file="bglog_12_22.RData",compress=TRUE) 
