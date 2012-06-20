rm(list=ls())
source("pphyproevo.R")

## test run from 1:10
system.time(res <- MLE.s(c(beta[1],gamma[1]),noNAind))
save(res,file="test2")

