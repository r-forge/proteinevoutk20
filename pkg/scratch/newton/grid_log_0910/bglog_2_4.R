rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4.8,-6.52),1:106)) 
 save(res,file="bglog_2_4.RData",compress=TRUE) 
