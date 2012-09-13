rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.8,-3.32),1:106)) 
 save(res,file="bglog_2_24.RData",compress=TRUE) 
