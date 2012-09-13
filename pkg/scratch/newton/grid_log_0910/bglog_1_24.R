rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-3.32),1:106)) 
 save(res,file="bglog_1_24.RData",compress=TRUE) 
