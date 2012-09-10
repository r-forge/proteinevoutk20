rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-12),1:106)) 
 save(res,file="bglog_12_3.RData",compress=TRUE) 
