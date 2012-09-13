rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4,-4.12),1:106)) 
 save(res,file="bglog_6_19.RData",compress=TRUE) 
