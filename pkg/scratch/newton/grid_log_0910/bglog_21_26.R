rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-3),1:106)) 
 save(res,file="bglog_21_26.RData",compress=TRUE) 
