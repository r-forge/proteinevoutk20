rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.6,-3),1:106)) 
 save(res,file="bglog_18_26.RData",compress=TRUE) 
