rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-1.6,-5.56),1:106)) 
 save(res,file="bglog_18_10.RData",compress=TRUE) 
