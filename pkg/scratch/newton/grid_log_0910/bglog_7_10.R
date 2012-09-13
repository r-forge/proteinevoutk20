rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-3.8,-5.56),1:106)) 
 save(res,file="bglog_7_10.RData",compress=TRUE) 
