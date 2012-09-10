rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-12,-13),1:106)) 
 save(res,file="bglog_3_2.RData",compress=TRUE) 
