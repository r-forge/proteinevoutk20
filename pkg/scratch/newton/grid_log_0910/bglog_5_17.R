rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.2,-4.44),1:106)) 
 save(res,file="bglog_5_17.RData",compress=TRUE) 
