rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.399999999999999,-5.72),1:106)) 
 save(res,file="bglog_24_9.RData",compress=TRUE) 
