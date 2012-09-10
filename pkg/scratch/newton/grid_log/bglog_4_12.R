rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-11,-3),1:106)) 
 save(res,file="bglog_4_12.RData",compress=TRUE) 
