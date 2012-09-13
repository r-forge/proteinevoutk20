rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-5,-5.56),1:106)) 
 save(res,file="bglog_1_10.RData",compress=TRUE) 
