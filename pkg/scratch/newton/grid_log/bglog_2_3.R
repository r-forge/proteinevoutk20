rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-13,-12),1:106)) 
 save(res,file="bglog_2_3.RData",compress=TRUE) 
