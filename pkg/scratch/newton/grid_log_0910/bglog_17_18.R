rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.8,-4.28),1:106)) 
 save(res,file="bglog_17_18.RData",compress=TRUE) 
