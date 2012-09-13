rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.4,-5.4),1:106)) 
 save(res,file="bglog_4_11.RData",compress=TRUE) 
