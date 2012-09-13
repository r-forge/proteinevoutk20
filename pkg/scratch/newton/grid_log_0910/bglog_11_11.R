rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-5.4),1:106)) 
 save(res,file="bglog_11_11.RData",compress=TRUE) 
