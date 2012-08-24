source("~/proteinevoutk20/pkg/R/pphyproevo.R")
load("~/proteinevoutk20/pkg/scratch/newton/sGrid.RData")
system.time(mle <- MLE.bg(c(0.1,0.01),c(0,0),c(1,1),s.lt1.index,trace=1,multicore=TRUE))

save(mle,file="mle.slt1.RData")
