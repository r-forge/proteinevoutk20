rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-12,-14),1:106)) 
 save(res,file="bglog_3_1.RData",compress=TRUE) 
