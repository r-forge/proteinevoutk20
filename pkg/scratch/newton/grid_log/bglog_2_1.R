rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-13,-14),1:106)) 
 save(res,file="bglog_2_1.RData",compress=TRUE) 
