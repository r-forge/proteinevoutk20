rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-10,-14),1:106)) 
 save(res,file="bglog_5_1.RData",compress=TRUE) 
