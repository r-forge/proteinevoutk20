rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-3.48),1:106)) 
 save(res,file="bglog_1_23.RData",compress=TRUE) 
