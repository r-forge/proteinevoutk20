rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.2,-5.56),1:106)) 
 save(res,file="bglog_15_10.RData",compress=TRUE) 
