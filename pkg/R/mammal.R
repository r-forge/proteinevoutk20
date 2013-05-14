source("~/proteinevoutk20/pkg/R/prune.R")
gene = 1
fastafile <- "~/proteinevoutk20/pkg/Data/mammals/mam15.fasta"
load("~/proteinevoutk20/pkg/Data/mammals/charset.RData") 
mamtree <- read.nexus("proteinevoutk20/pkg/Data/mammals/T15taxa.nex")
prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/mammal/mam15_",gene,"_prottest.txt",sep="")
RDatafile <- paste("gene",gene,".RData",sep="")
best_emp_model <- get_best_model(prottestfile)
dtip = 12
p2 <- prune_emp(fastafile,dtip,mamtree,best_emp_model$model,range=charset[[gene]])
p1 <- prune_new(fastafile,dtip,mamtree,ancestral="max",range=charset[[gene]])
save.image(RDatafile,compress=TRUE)
##truncate the names of tip so that the length is no bigger than 10 for prottest
## the resulting tree is written in file mam15.tre
## mamtree$tip.label <- sapply(1:15, function(x) substr(mamtree$tip.label[x],1,10))


# for(i in 1:97){
#   mam <- conv("~/proteinevoutk20/pkg/Data/mammals/mam15.fasta",range=charset[[i]],"AA")
#   mamlist <- lapply(seq_len(nrow(mam)),function(i) mam[i,])
#   names(mamlist) <- dimnames(mam)[[1]]
#   filename = paste("mam_15_",i,".nex",sep="")
#   write.nexus.data(mamlist,file=filename,format="protein",interleaved=FALSE)
# }