## beetles data analysis
## partitions: 
## AATS: 1-915
## CAD: 916-2962
## TPI: 2963-3460 
## PGD: 3461-4264 
## SNF: 4265-4825 
## RNA POL II: 4826-5725
## EF-1alpha-1,2,3pos: 5726-6783
## 18S: 6784-10478
## 28S: 10479-12777
## for data without 18S and 28S
pos1 <- 1:915
pos2 <- 916:2973
pos3 <- 2974:3472
pos4 <- 3473:4275
pos5 <- 4276:4836
pos6 <- 4837:5736
bedata=beetle
for(i in 1:length(beetle))
  bedata[[i]] = beetle[[i]][pos4]
bedata <- conv(beetle,range=pos6,type="AA")
sum(bedata=="*")
source("~/proteinevoutk20/pkg/R/main.R")
#beetle data without 18S and 28S, 28 taxa from earlier paper
beetle <- read.nexus.data("~/proteinevoutk20/pkg/Data/beetle_short_data.nex")
#beetle tree, 28 taxa
betree <- read.nexus("~/proteinevoutk20/pkg/Data/beetle_short.nex")
betree$edge.length <- rep(0.001,56)
pos <- 1:915
bedata <- conv(beetle,range=pos,type="phyDat")
be.res <- mllm1(bedata,betree,s=1,beta=be,gamma=ga)
be.optim <- optim.mllm1(be.res,optQ=T,optBranch=T,optsWeight=T,
                        control=list(epsilon=1e-08,hmaxit=30,htrace=TRUE,print_level=0,maxit=100))
save(be.optim,file="beetle_optim.RData")
#parts <- c(0,915,2973,3471,4275,4836,5736,6783,10478,12777)
#partlens <- (parts[-1]-parts[-10])/3
