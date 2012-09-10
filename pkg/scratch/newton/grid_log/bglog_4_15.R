rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-11,0),1:106)) 
 save(res,file="bglog_4_15.RData",compress=TRUE) 
