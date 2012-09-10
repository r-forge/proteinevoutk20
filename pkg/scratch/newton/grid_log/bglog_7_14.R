rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-8,-1),1:106)) 
 save(res,file="bglog_7_14.RData",compress=TRUE) 
