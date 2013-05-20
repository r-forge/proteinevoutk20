#char <- list(aats=1:915,cad=916:2962,tpi=2963:3460,pgd=3461:4264,
 #               snf=4265:4825,pol=4826:5725,ef1a=5726:6783,whole=1:6783)

test = TRUE
beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetle/beetles.fasta")
load("~/proteinevoutk20/pkg/Data/beetle/charset.RData") 
beetree <- read.nexus("~/proteinevoutk20/pkg/Data/beetle/tree1.nex")
fastafile = beetle
gene = 6
prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/beetle/beetle34_",gene,"_prottest.txt",sep="")
RDatafile <- paste("gene",gene,".RData",sep="")
best_emp_model <- get_best_model(prottestfile)
#dtip = "Trichoptera.Hydropsyche"
dtip = "Mecoptera.Nannochorista"
p2 <- prune_emp(fastafile,dtip,beetree,best_emp_model$model,range=char[[gene]])
p1 <- prune_new(fastafile,dtip,beetree,ancestral="max",range=char[[gene]],ancestral="eqm")
save.image(RDatafile,compress=TRUE)

