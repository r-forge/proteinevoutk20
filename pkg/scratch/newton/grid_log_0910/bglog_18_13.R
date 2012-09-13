rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.6,-5.08),1:106)) 
 save(res,file="bglog_18_13.RData",compress=TRUE) 
