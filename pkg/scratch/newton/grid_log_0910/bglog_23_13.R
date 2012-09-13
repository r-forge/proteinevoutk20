rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.6,-5.08),1:106)) 
 save(res,file="bglog_23_13.RData",compress=TRUE) 
