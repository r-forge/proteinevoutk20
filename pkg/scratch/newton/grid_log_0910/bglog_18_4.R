rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.6,-6.52),1:106)) 
 save(res,file="bglog_18_4.RData",compress=TRUE) 
