rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-4.28),1:106)) 
 save(res,file="bglog_21_18.RData",compress=TRUE) 
