rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-4.28),1:106)) 
 save(res,file="bglog_26_18.RData",compress=TRUE) 
