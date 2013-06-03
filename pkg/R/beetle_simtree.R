source("~/proteinevoutk20/pkg/R/prune.R") 

gene = 6
beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetle/beetles.fasta")
load("~/proteinevoutk20/pkg/Data/beetle/charset.RData") 
beetree <- read.nexus("~/proteinevoutk20/pkg/Data/beetle/tree1.nex")
fastafile = beetle
prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/beetle/beetle34_",gene,"_prottest.txt",sep="")
best_emp_model <- get_best_model(prottestfile)
RDatafile <- paste("~/BackupProEvo/Newton/beetle/tree1/rootMax/gene",gene,".RData",sep="")
imagefile <- paste("gene",gene,".RData",sep="") 
best_emp_model <- get_best_model(prottestfile) 
data = conv(fastafile,type="phyDat")
load(RDatafile)
#######################################################################
## simulation on one branch only
# test = FALSE 
# p2 <- prune_emp(fastafile,"Smik",ROKAS_TREE,best_emp_model$model) 
# p1 <- prune_new(fastafile,"Smik",ROKAS_TREE,ancestral="eqm") 
# save.image(RDatafile,compress=TRUE) 
#######################################################################
## simulation along the tree
res <- mllm1(res_op$data,res_op$tree,Qall=res_op$Qall,ancestral="max")
nr <- attr(data,"nr")
index <- attr(data,"index")
opaa <- res$ll$opaa
s = res_op$s
beta = res_op$GMweights[2]
gamma = res_op$GMweights[3]
nsim <- 10
sim <- vector(mode="list")
sim$new <- vector(mode="list",length=nsim)
sim$emp <- vector(mode="list",length=nsim)
for(i in 1:nsim){
  if(exists(".Random.seed"))
    rm(.Random.seed)
  ## root sequence, from the inference of ancestral states in mllm1, length equal to number of sites
  rootseq <- sapply(1:length(index),function(x) sample(1:20,size=1,replace=TRUE,prob=res$ll$ancestral[,index[x]]))
  sim$new[[i]] <- simTree(tree=res$tree,protein_op=opaa[index],matall=res$Qall,rootseq=rootseq)
  sim$emp[[i]] <- sim_emp(data,best_emp_model) ## best empirical model
}
save.image(file=imagefile,compress=TRUE)