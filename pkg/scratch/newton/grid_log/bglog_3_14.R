rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-12,-1),1:106)) 
 save(res,file="bglog_3_14.RData",compress=TRUE) 
