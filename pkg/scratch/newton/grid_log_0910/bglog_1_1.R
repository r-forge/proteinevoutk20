rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-7),1:106)) 
 save(res,file="bglog_1_1.RData",compress=TRUE) 
