rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-8,-11),1:106)) 
 save(res,file="bglog_7_4.RData",compress=TRUE) 
