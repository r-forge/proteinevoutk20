rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-14,-3),1:106)) 
 save(res,file="bglog_1_12.RData",compress=TRUE) 
