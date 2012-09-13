rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.4,-4.28),1:106)) 
 save(res,file="bglog_9_18.RData",compress=TRUE) 
