rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-4.6),1:106)) 
 save(res,file="bglog_26_16.RData",compress=TRUE) 
