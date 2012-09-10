rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3,-13),1:106)) 
 save(res,file="bglog_12_2.RData",compress=TRUE) 
