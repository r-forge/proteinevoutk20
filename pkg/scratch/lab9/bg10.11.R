rm(list=ls()) 
 source("pphyproevo.R") 
 system.time(res <- MLE.s(c(beta[10],gamma[11]),1:106)) 
 save(res,file="bg.10.11",compress=TRUE) 
