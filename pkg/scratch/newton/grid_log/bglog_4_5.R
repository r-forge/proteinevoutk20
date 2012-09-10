rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-11,-10),1:106)) 
 save(res,file="bglog_4_5.RData",compress=TRUE) 
