rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-9,-2),1:106)) 
 save(res,file="bglog_6_13.RData",compress=TRUE) 
