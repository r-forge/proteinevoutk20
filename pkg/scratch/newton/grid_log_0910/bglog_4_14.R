rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.4,-4.92),1:106)) 
 save(res,file="bglog_4_14.RData",compress=TRUE) 
