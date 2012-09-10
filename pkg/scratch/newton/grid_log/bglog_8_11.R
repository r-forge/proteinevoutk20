rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-7,-4),1:106)) 
 save(res,file="bglog_8_11.RData",compress=TRUE) 
