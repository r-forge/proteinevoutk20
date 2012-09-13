rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-7),1:106)) 
 save(res,file="bglog_11_1.RData",compress=TRUE) 
