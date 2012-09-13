rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.8,-5.88),1:106)) 
 save(res,file="bglog_22_8.RData",compress=TRUE) 
