rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-14,-2),1:106)) 
 save(res,file="bglog_1_13.RData",compress=TRUE) 
