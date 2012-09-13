rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.6,-6.52),1:106)) 
 save(res,file="bglog_3_4.RData",compress=TRUE) 
