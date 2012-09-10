rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-11,-7),1:106)) 
 save(res,file="bglog_4_8.RData",compress=TRUE) 
