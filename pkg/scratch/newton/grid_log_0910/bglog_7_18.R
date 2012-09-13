rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.8,-4.28),1:106)) 
 save(res,file="bglog_7_18.RData",compress=TRUE) 
