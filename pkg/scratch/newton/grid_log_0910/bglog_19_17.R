rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.4,-4.44),1:106)) 
 save(res,file="bglog_19_17.RData",compress=TRUE) 
