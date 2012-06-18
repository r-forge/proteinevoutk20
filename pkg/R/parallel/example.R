source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
data(geneData,package="Biobase")
fakeData <- cbind(geneData,geneData,geneData,geneData)
library(boot)
pair <- combn(1:nrow(geneData),2,simplify=F)

geneCor2 <- function(x, gene = fakeData) {
  mydata <- cbind(gene[x[1], ], gene[x[2], ])
  mycor <- function(x, i) cor(x[i,1], x[i,2])
  boot.out <- boot(mydata, mycor, 1000)
  boot.ci(boot.out, type = "bca")$bca[4:5]
}
system.time(out <- lapply(pair[1:10],geneCor2))
pair2 <- sample(pair,300)
system.time(out <- lapply(pair2,geneCor2))
system.time(out <- mclapply(pair2,geneCor2))

library(snow)
hosts <- c(rep("localhost",2),rep("sisal",8))
cl <- makeCluster(hosts,type="SOCK")
clusterCall(cl,date)
clusterCall(cl,function(x) Sys.info()[c("nodename","machine")])
clusterExport(cl,"fakeData")
clusterEvalQ(cl,library(boot))
system.time(out <- clusterApplyLB(cl,pair2,geneCor2))
tm <- snow.time(out <- clusterApplyLB(cl,pair2,geneCor2))
stopCluster(cl)