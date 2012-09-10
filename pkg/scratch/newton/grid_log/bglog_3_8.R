rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-12,-7),1:106)) 
 save(res,file="bglog_3_8.RData",compress=TRUE) 
