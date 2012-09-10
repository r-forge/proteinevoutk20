rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-13,-4),1:106)) 
 save(res,file="bglog_2_11.RData",compress=TRUE) 
