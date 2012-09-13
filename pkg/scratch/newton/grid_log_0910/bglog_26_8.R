rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(0,-5.88),1:106)) 
 save(res,file="bglog_26_8.RData",compress=TRUE) 
