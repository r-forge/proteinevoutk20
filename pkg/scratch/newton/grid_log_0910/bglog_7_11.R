rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.8,-5.4),1:106)) 
 save(res,file="bglog_7_11.RData",compress=TRUE) 
