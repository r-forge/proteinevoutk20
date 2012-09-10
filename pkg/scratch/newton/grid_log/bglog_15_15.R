rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,0),1:106)) 
 save(res,file="bglog_15_15.RData",compress=TRUE) 
