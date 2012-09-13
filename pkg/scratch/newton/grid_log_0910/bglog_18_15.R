rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.6,-4.76),1:106)) 
 save(res,file="bglog_18_15.RData",compress=TRUE) 
