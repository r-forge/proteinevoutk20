#tvec <- c("010","025","050","1") #parameters to pass into the function
tvec <- c("010")
line1 <- "rm(list=ls())"
line2 <- 'source("simsumPHYLOGENfunc.R")'
pbsCommands=""
queue="short*"

for(j in sequence(length(tvec))){
  line3 <- paste('x <- read.tree("trees100.',tvec[j],'")',sep="")
  for(i in 1:1){
    lsString=paste(paste('ls -1 testB*',' | grep -c testB',i,'.100.',tvec[j],sep="",collapse=""))
    matchCount=suppressWarnings(as.numeric(system(lsString,intern=TRUE)))
    if (matchCount!=1) {
      line4 <- paste("test.",i,".100.",tvec[j]," <- Null.null(x,start=",i,",stop=",i,",reps=1000,mod.id=c(1,1,1,1))",sep="")
      line5 <- paste('save(test.',i,'.100.',tvec[j],',file="testB',i,'.100.',tvec[j],'",compress=TRUE)',sep="")
      cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",file=paste("CMD100_B",i,".",tvec[j],".R",sep=""))
      pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
      pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
      pbsCommands=paste(pbsCommands,'#$ -M omeara.brian@gmail.com', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
      pbsCommands=paste(pbsCommands,"\n","#$ -N B",i,"J",tvec[j],sep="")
      pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH CMD100_B",i,".",tvec[j],".R",sep="")
      cat(pbsCommands,file=paste('runB',i,'J',tvec[j],'.sh',sep=""),append=FALSE)
      system(paste('chmod u+x runB',i,'J',tvec[j],'.sh',sep=""))
      system(paste('/opt/sge/bin/lx24-amd64/qsub runB',i,'J',tvec[j],'.sh',sep=""))
    }
     while(as.numeric(system("/opt/sge/bin/lx24-amd64/qstat | grep -c bomeara",intern=TRUE))>1000) {
       Sys.sleep(37)
    }
  }
}

