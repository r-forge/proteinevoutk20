rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-2),1:106)) 
 save(res,file="bglog_15_13.RData",compress=TRUE) 
