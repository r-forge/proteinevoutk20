##parameter (beta, gamma) values to run
# l <- 100
# beta <- seq(0,0.5,length.out=(l+1))[-1]
# gamma <- seq(0,0.05,length.out=(l+1))[-1]
##rm(list=ls())
line1 <- "rm(list=ls())"
##source("pphyproevo.R")
line2 <- 'source("pphyproevo.R")'
pbsCommands=""
queue="short*"

for(i in 1:l){
  for(j in 1:l){
    for(gene in 1:106){
      ##gdata <- data[[gene]]
      line3 <- paste("gdata <- data[[",gene,"]]",sep="")
      ##system.time(res <- MLE_GTR(1,0,1e4,tree,gdata,al,beta[i],gamma[j],mumat,optim.m=1))
      line4 <- paste("system.time(res <- MLE_GTR(1,0,1e5,tree,gdata,al,beta[",i,"],gamma[",j,"],mumat,optim.m=1))",sep="")
      ##save(res,file="weight.i.j.gene",compress=TRUE)
      line5 <- paste('save(res,file="gene.',i,'.',j,'.',gene,'",compress=TRUE)',sep="")
      cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",file=paste("gene",i,".",j,".",gene,".R",sep=""))
      pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
      pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
      pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
      pbsCommands=paste(pbsCommands,"\n","#$ -N B",i,"_",j,"_",gene,sep="")
      pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH gene",i,".",j,".",gene,".R",sep="")
      cat(pbsCommands,file=paste('runB',i,'_',j,'_',gene,'.sh',sep=""),append=FALSE)
      system(paste('chmod u+x runB',i,'_',j,'_',gene,'.sh',sep=""))
      system(paste('/opt/sge/bin/lx24-amd64/qsub runB',i,'_',j,'_',gene,'.sh',sep=""))
    }
    while(as.numeric(system("/opt/sge/bin/lx24-amd64/qstat | grep -c jchai1",intern=TRUE))>1000) {
      Sys.sleep(37)
  }
}
