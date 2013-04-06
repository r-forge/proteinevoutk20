charset <- list(aats=1:915,cad=916:2962,tpi=2963:3460,pgd=3461:4264,
                snf=4265:4825,pol=4826:5725,ef1a=5726:6783,whole=1:6783)
beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetles.fasta")
bedata <- conv(beetle,range=charset[[7]],type="phyDat",frame=2)
sort(unique(c(bedata)))
sum(bedata=="*")

ch = charset[[7]]
sapply(1:34,function(i) sum(beetle[[i]][ch]=="-"))
cat(beetle[[2]][ch],sep="")

data = beetle[[7]][ch]
data = data[which(data!="-")]

bdata <- conv(beetle,range=char[[4]],type="AA",frame=0)
sort(unique(c(bdata)))
sum(bdata=="*")

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


source("~/proteinevoutk20/pkg/R/main.R") 
load("~/proteinevoutk20/pkg/Data/beetle/charset.RData") 
beetree <- read.nexus("~/proteinevoutk20/pkg/Data/beetle/tree1.nex") 
range = char[[1]] 
beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetle/beetles.fasta") 
bee = conv(beetle, range=range,"phyDat") 
sort(unique(unlist(bee))) 
bee.res = mllm1(bee,beetree,s=1,beta=be,gamma=ga,Q=rep(1,6),Ne=700) 
bee.res$ll$loglik 
bee.optim = optim.mllm(bee.res,optQ=T,optBranch=T,optsWeight=T,optOpw=FALSE,control=list(epsilon=1e-08,hmaxit=1,htrace=1,print_level=0,maxeval="20"),Ne=700) 
save(bee.optim,file="bee_1.RData",compress=TRUE) 