rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-9,-7),1:106)) 
 save(res,file="bglog_6_8.RData",compress=TRUE) 
