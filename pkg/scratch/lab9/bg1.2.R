rm(list=ls()) 
 source("~/proteinevoutk20/pkg/R/pphyproevo.R") 
 system.time(res <- MLE.s(c(beta[1],gamma[2]),1:106)) 
 save(res,file="bg.1.2",compress=TRUE) 
