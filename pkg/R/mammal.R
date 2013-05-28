source("~/proteinevoutk20/pkg/R/prune.R")
gene = 1
fastafile <- "~/proteinevoutk20/pkg/Data/mammals/mam15.fasta"
load("~/proteinevoutk20/pkg/Data/mammals/charset.RData") 
mamtree <- read.nexus("~/proteinevoutk20/pkg/Data/mammals/T15taxa.nex")
prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/mammal/mam15_",gene,"_prottest.txt",sep="")
RDatafile <- paste("gene",gene,".RData",sep="")
best_emp_model <- get_best_model(prottestfile)
data = conv(fastafile,range=charset[[gene]],type="phyDat")
#######################################################################
## simulation along one branch
dtip = 8 #"Ceratotherium_simum"
p2 <- prune_emp(fastafile,dtip,mamtree,best_emp_model$model,range=charset[[gene]])
p1 <- prune_new(fastafile,dtip,mamtree,ancestral="eqm",range=charset[[gene]])
#######################################################################
## simulation along the tree
nsites = length(charset[[gene]])/3
model <- best_emp_model$model
is.G <- "G" %in% model #gamma?
is.F <- "F" %in% model #empirical frequencies?
is.I <- "I" %in% model #invariant sites included?
bf = NULL
k = 1
shape = 1
if(is.F) bf=findBf2(data)
if(is.I) inv=best_emp_model$inv
if(is.G) {shape = best_emp_model$shape 
          k=4}
sim_emp <- simTreeEmp(best_emp_model$tree,l=nsites,bf=bf,rate=shape,k=k,model=model[1])