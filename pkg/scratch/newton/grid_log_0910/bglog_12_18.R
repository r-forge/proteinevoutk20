rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.8,-4.28),1:106)) 
 save(res,file="bglog_12_18.RData",compress=TRUE) 
