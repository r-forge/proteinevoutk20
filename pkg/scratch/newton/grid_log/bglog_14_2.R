rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1,-13),1:106)) 
 save(res,file="bglog_14_2.RData",compress=TRUE) 
