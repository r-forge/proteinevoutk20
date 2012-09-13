rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-4.76),1:106)) 
 save(res,file="bglog_11_15.RData",compress=TRUE) 
