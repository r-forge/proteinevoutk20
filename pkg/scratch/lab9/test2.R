rm(list=ls())
source("~/proteinevoutk20/pkg/R/pphyproevo.R")

## test run from 11:106
system.time(res <- MLE.s(c(beta[1],gamma[1]),11:106))
save.image(file="test2")
