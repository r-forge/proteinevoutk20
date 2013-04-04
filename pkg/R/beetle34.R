# charset aats: 1-915;
# charset cad: 916-2962;
# charset tpi: 2963-3460;
# charset pgd:  3461-4264;
# charset snf: 4265-4825;
# charset pol: 4826-5725;
# charset ef1a: 5726-6783;
# charset 18s: 6784-10478;
# charset 28s: 10479-12777;
charset <- list(aats=1:915,cad=916:2962,tpi=2963:3460,pgd=3461:4264,
                snf=4265:4825,pol=4826:5725,ef1a=5726:6783,whole=1:6783)
beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetles.fasta")
bedata <- conv(beetle,range=charset[[2]],type="AA",frame=2)
sort(unique(c(bedata)))
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