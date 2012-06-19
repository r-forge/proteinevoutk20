rm(list=ls()) 
 source("~/proteinevoutk20/pkg/R/pphyproevo.R") 
 system.time(res <- MLE.s(c(beta[3],gamma[2]),1:106)) 
 save(res,file="bg.3.2",compress=TRUE) 
