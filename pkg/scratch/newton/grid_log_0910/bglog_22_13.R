rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.8,-5.08),1:106)) 
 save(res,file="bglog_22_13.RData",compress=TRUE) 
