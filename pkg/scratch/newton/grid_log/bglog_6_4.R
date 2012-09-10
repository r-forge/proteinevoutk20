rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-9,-11),1:106)) 
 save(res,file="bglog_6_4.RData",compress=TRUE) 
