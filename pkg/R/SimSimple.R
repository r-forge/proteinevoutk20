sites <- c(100,300,500,700,1000)
s <- c(0.1,0.3,0.5,0.7,0.9,1)
tip <- c(4,4,8,8,12,12)

line1 <- 'load("~/proteinevoutk20/pkg/R/simtrees2.RData")'
line1.1 <- 'source("~/proteinevoutk20/pkg/R/simulation.R")'
line1.2 <- 'source("~/proteinevoutk20/pkg/R/hessian.R")'
pbsCommands=""
queue="long*"

for(k in 1:6){
  for(i in sites){
    for(j in s){
      line2 <- paste("start_seq=sample(20,",i,",replace=T)",sep="")
      line3 <- paste("op_seq=sample(20,",i,",replace=T)",sep="")
      line4 <- paste("system.time(sim <- simTree(trees[[",k,"]]",",",i,",op_seq,20,",j,",Nu_vec,alpha=1,beta=0.5,gamma=0.1,rootseq=start_seq, ancestral=TRUE, simple=TRUE))",sep="")
      line5 <- paste("system.time(mle <- MLE_sw(c(0.1,0.1,0.1),rep(0,3),rep(1,3),trees[[",k,"]],sim[1:",tip[k],",],m=20,1,mumat,protein_op=op_seq,root=start_seq,trace=1))",sep="")
      line6 <- paste("hes <- find_hessian(mle$par,trees[[",k,"]],sim[1:tip[",k,"],],1,mumat,20,protein_op=op_seq,root=start_seq)",sep="")
      line7 <- paste('save.image(file="~/proteinevoutk20/pkg/RData/sim_simple/simSimple_',k,'_',i,'_',j,'.RData",compress=TRUE)',sep="")
      
      cat(line1,"\n",line1.1,"\n",line1.2,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",line6, "\n",line7, "\n",file=paste("simSimple_",k,"_",i,"_",j,".R",sep=""))
      ##bashCommands=paste(bashCommands,"nohup R CMD BATCH --vanilla sim_",k,"_",i,"_",j,".R > /dev/null &","\n",sep="")
      ##cat(bashCommands, file="simjobssub.sh")
      pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
      pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
      pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
      pbsCommands=paste(pbsCommands,"\n","#$ -N SimSimple_",k,"_",i,"_",j,sep="")
      pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH simSimple_",k,"_",i,"_",j,".R",sep="")
      cat(pbsCommands,file=paste('runSimSimple_',k,'_',i,'_',j,'.sh',sep=""),append=FALSE)
      system(paste('chmod u+x runSimSimple_',k,'_',i,'_',j,'.sh',sep=""))
      system(paste('/opt/sge/bin/lx24-amd64/qsub runSimSimple_',k,'_',i,'_',j,'.sh',sep=""))
    }
  }
}
