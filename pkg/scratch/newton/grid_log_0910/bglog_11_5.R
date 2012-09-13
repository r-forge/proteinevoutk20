rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-6.36),1:106)) 
 save(res,file="bglog_11_5.RData",compress=TRUE) 
