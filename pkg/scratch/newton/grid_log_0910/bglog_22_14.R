rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.8,-4.92),1:106)) 
 save(res,file="bglog_22_14.RData",compress=TRUE) 
