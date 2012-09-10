rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-8,-9),1:106)) 
 save(res,file="bglog_7_6.RData",compress=TRUE) 
