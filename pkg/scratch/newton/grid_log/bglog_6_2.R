rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-9,-13),1:106)) 
 save(res,file="bglog_6_2.RData",compress=TRUE) 
