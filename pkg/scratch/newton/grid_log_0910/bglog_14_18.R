rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.4,-4.28),1:106)) 
 save(res,file="bglog_14_18.RData",compress=TRUE) 
