rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-9,-3),1:106)) 
 save(res,file="bglog_6_12.RData",compress=TRUE) 
