rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.6,-5.56),1:106)) 
 save(res,file="bglog_13_10.RData",compress=TRUE) 
