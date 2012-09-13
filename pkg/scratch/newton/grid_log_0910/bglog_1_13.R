rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-5.08),1:106)) 
 save(res,file="bglog_1_13.RData",compress=TRUE) 
