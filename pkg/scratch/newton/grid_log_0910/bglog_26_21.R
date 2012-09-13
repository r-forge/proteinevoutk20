rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-3.8),1:106)) 
 save(res,file="bglog_26_21.RData",compress=TRUE) 
