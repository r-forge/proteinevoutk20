rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.2,-5.56),1:106)) 
 save(res,file="bglog_10_10.RData",compress=TRUE) 
