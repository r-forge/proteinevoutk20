rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4,-1),1:106)) 
 save(res,file="bglog_11_14.RData",compress=TRUE) 
