rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-4.92),1:106)) 
 save(res,file="bglog_1_14.RData",compress=TRUE) 
