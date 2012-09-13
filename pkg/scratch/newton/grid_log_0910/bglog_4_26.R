rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.4,-3),1:106)) 
 save(res,file="bglog_4_26.RData",compress=TRUE) 
