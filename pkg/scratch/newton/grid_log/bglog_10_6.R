rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-9),1:106)) 
 save(res,file="bglog_10_6.RData",compress=TRUE) 
