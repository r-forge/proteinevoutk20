rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-12,-5),1:106)) 
 save(res,file="bglog_3_10.RData",compress=TRUE) 
