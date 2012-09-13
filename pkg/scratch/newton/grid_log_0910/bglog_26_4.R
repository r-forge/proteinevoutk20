rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-6.52),1:106)) 
 save(res,file="bglog_26_4.RData",compress=TRUE) 
