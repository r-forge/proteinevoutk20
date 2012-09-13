rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s_log(c(-0.399999999999999,-4.12),1:106)) 
 save(res,file="bglog_24_19.RData",compress=TRUE) 
