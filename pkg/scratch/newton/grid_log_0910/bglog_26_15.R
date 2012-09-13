rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-4.76),1:106)) 
 save(res,file="bglog_26_15.RData",compress=TRUE) 
