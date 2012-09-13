rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.4,-5.88),1:106)) 
 save(res,file="bglog_14_8.RData",compress=TRUE) 
