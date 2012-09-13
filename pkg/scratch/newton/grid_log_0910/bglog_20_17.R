rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.2,-4.44),1:106)) 
 save(res,file="bglog_20_17.RData",compress=TRUE) 
