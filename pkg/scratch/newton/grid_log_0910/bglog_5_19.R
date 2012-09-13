rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-4.12),1:106)) 
 save(res,file="bglog_5_19.RData",compress=TRUE) 
