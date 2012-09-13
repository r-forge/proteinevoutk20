rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.8,-3.16),1:106)) 
 save(res,file="bglog_22_25.RData",compress=TRUE) 
