rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.4,-6.68),1:106)) 
 save(res,file="bglog_4_3.RData",compress=TRUE) 
