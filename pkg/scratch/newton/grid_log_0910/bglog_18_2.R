rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.6,-6.84),1:106)) 
 save(res,file="bglog_18_2.RData",compress=TRUE) 
