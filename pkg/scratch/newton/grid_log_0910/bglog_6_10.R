rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-4,-5.56),1:106)) 
 save(res,file="bglog_6_10.RData",compress=TRUE) 
