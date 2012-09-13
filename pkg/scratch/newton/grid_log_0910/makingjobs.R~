##############################################
##parameter (beta, gamma) values to run
#l <- 20
# beta <- seq(0,1,length.out=(l+1))[-1]
# gamma <- seq(0,1,length.out=(l+1))[-1]

l <- 15
beta <- seq(-15,0,length.out=(l+1))[-1]
gamma <- seq(-15,0,length.out=(l+1))[-1]
##############################################
##rm(list=ls())
line1 <- "rm(list=ls())"
##source("pphyproevo.R")
line2 <- 'source("pphyproevo.R")'
pbsCommands=""
queue="medium*"

for(i in 1:l){
  for(j in 1:l){
      ##system.time(res <- MLE.s(c(beta[i],gamma[j]),1:106))
      line3 <- paste("system.time(res <- MLE.s_log(c(",beta[i],",",gamma[j],"),1:106))",sep="")
      ##save(res,file="bg.i.j",compress=TRUE)
      line4 <- paste('save(res,file="bglog_',i,'_',j,'.RData",compress=TRUE)',sep="")
      cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",file=paste("bglog_",i,"_",j,".R",sep=""))
      pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
      pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
      pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
      pbsCommands=paste(pbsCommands,"\n","#$ -N B_",i,"_",j,sep="")
      pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH bglog_",i,"_",j,".R",sep="")
      cat(pbsCommands,file=paste('runB_',i,'_',j,'.sh',sep=""),append=FALSE)
      system(paste('chmod u+x runB_',i,'_',j,'.sh',sep=""))
      system(paste('/opt/sge/bin/lx24-amd64/qsub runB_',i,'_',j,'.sh',sep=""))
    }
     while(as.numeric(system("/opt/sge/bin/lx24-amd64/qstat | grep -c jchai1",intern=TRUE))>1000) {
       Sys.sleep(37)
   }
}
