rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-5.4),1:106)) 
 save(res,file="bglog_1_11.RData",compress=TRUE) 
