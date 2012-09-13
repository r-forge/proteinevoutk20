rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.6,-6.68),1:106)) 
 save(res,file="bglog_23_3.RData",compress=TRUE) 
