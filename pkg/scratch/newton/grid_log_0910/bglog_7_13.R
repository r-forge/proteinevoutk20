rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.8,-5.08),1:106)) 
 save(res,file="bglog_7_13.RData",compress=TRUE) 
