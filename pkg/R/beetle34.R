#charset <- list(aats=1:915,cad=916:2962,tpi=2963:3460,pgd=3461:4264,
 #               snf=4265:4825,pol=4826:5725,ef1a=5726:6783,whole=1:6783)

test = TRUE

beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetle/beetles.fasta")
load("~/proteinevoutk20/pkg/Data/beetle/charset.RData") 
beetree <- read.nexus("~/proteinevoutk20/pkg/Data/beetle/tree1.nex")
fastafile = beetle

gene = 1
prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/beetle/beetle34_",gene,"_prottest.txt",sep="")
RDatafile <- paste("gene",gene,".RData",sep="")
best_emp_model <- get_best_model(prottestfile)
dtip = 17
p2 <- prune_emp(fastafile,dtip,beetree,best_emp_model$model,range=char[[gene]])
p1 <- prune_new(fastafile,dtip,beetree,ancestral="max",range=char[[gene]])
save.image(RDatafile,compress=TRUE)


beetle <- seqinr::read.fasta("~/proteinevoutk20/pkg/Data/beetle/beetles.fasta")
load("~/proteinevoutk20/pkg/Data/beetle/charset.RData") 
for(i in 1:8){
  bee <- conv(beetle,range=char[[i]],"AA")
  beelist <- lapply(seq_len(nrow(bee)),function(i) bee[i,]) #turn matrix into list
  names(beelist) <- beename
  filename = paste("beetle34_",i,"_shortname.nex",sep="")
  write.nexus.data(beelist,file=filename,format="protein",interleaved=FALSE)
}

ch = charset[[7]]
sapply(1:34,function(i) sum(beetle[[i]][ch]=="-"))
cat(beetle[[2]][ch],sep="")

data = beetle[[7]][ch]
data = data[which(data!="-")]

bdata <- conv(beetle,range=char[[4]],type="AA",frame=0)
sort(unique(c(bdata)))
sum(bdata=="*")
