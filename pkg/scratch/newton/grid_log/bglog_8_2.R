rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-7,-13),1:106)) 
 save(res,file="bglog_8_2.RData",compress=TRUE) 
