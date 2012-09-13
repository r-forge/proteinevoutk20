rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-5.24),1:106)) 
 save(res,file="bglog_11_12.RData",compress=TRUE) 
