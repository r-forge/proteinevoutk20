rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-3.48),1:106)) 
 save(res,file="bglog_11_23.RData",compress=TRUE) 
