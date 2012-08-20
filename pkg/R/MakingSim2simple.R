##source("~/proteinevoutk20/pkg/R/simulation.R")


sites <- c(100,500,1000)
s <- c(0,0.01,0.1)

line1 <- 'source("~/proteinevoutk20/pkg/R/simulation.R")'
line1.1 <- 'load("~/proteinevoutk20/pkg/R/simtrees2.RData")'
pbsCommands=""
queue="medium*"

for(k in 1:6){
  for(i in sites){
    for(j in s){
      line2 <- paste("start_seq=sample(20,",i,",replace=T)",sep="")
      line3 <- paste("op_seq=rep(1,",i,")",sep="")
      line4 <- paste("system.time(sim <- simTree(trees[[",k,"]]",",",i,",op_seq,20,",j,",Nu_vec,rootseq=start_seq, ancestral=TRUE, simple=TRUE))",sep="")
      line5 <- paste('save.image(file="sim2simple_',k,'_',i,'_',j,'.RData",compress=TRUE)',sep="")
      cat(line1,"\n",line1.1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",file=paste("sim2simple_",k,"_",i,"_",j,".R",sep=""))
      ##bashCommands=paste(bashCommands,"nohup R CMD BATCH --vanilla sim_",k,"_",i,"_",j,".R > /dev/null &","\n",sep="")
      ##cat(bashCommands, file="simjobssub.sh")
      pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
      pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
      pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
      pbsCommands=paste(pbsCommands,"\n","#$ -N Sim2simple_",k,"_",i,"_",j,sep="")
      pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH sim2simple_",k,"_",i,"_",j,".R",sep="")
      cat(pbsCommands,file=paste('runSim2simple_',k,'_',i,'_',j,'.sh',sep=""),append=FALSE)
      system(paste('chmod u+x runSim2simple_',k,'_',i,'_',j,'.sh',sep=""))
      system(paste('/opt/sge/bin/lx24-amd64/qsub runSim2simple_',k,'_',i,'_',j,'.sh',sep=""))
    }
  }
}
