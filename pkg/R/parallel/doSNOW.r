#!/usr/bin/env Rscript
library(doSNOW)
hosts <- c( 'localhost','localhost', 'localhost','localhost','sisal','sisal')
cl <- makeCluster(hosts, type="SOCK")
registerDoSNOW(cl)
source("base.r")
stopCluster(cl)
