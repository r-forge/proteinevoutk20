rm(list=ls()) 
 source("~/proteinevoutk20/pkg/R/pphyproevo.R") 
 system.time(res <- MLE.s(c(beta[3],gamma[3]),1:106)) 
 save(res,file="bg.3.3",compress=TRUE) 
