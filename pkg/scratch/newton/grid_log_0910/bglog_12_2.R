rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-2.8,-6.84),1:106)) 
 save(res,file="bglog_12_2.RData",compress=TRUE) 
