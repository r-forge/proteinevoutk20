rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-4.92),1:106)) 
 save(res,file="bglog_21_14.RData",compress=TRUE) 
