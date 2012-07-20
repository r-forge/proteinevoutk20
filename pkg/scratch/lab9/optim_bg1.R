source("~/proteinevoutk20/pkg/R/pphyproevo.R")
system.time(mle <- MLE.bg(c(0.1,0.01),c(0,0),c(1,1),1:106,trace=1,multicore=TRUE))

save(mle,file="mlebg1.RData")
