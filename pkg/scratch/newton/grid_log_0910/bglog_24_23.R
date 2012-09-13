rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.399999999999999,-3.48),1:106)) 
 save(res,file="bglog_24_23.RData",compress=TRUE) 
