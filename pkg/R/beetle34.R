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
