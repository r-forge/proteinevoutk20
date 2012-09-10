rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-5),1:106)) 
 save(res,file="bglog_15_10.RData",compress=TRUE) 
