rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-5.24),1:106)) 
 save(res,file="bglog_26_12.RData",compress=TRUE) 
