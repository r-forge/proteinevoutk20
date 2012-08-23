##source("~/proteinevoutk20/pkg/R/simulation.R")


sites <- c(100,500,1000)
s <- c(0,0.01,0.1)
tip <- c(4,4,8,8,12,12)

line2 <- 'source("~/proteinevoutk20/pkg/R/pphyproevo.R")'
pbsCommands=""
queue="long*"
subjob=FALSE

for(k in 1:6){
  for(i in sites){
    for(j in s){
      ##load the workspace which stores the simulated data
      line1 <- paste('load("~/proteinevoutk20/pkg/RData/sim_notsimple/sim2_',k,'_',i,'_',j,'.RData")',sep="")
      ## find the mle for s, beta and gamma, given alpha value
      line3 <- paste("system.time(mle <- MLE_sw(c(0.1,0.01,0.001),rep(0,3),rep(1,3),trees[[",k,"]],sim[1:",tip[k],",],m=20,al,mumat,protein_op=op_seq,root=start_seq,trace=1))",sep="")
      line4 <- paste("hes <- find_hessian(mle[[1]]$par,trees[[",k,"]],sim[1:tip[",k,"],al,mumat,20,protein_op=op_seq,root=start_seq)",sep="")
      line5 <- paste('save.image(file="~/proteinevoutk20/pkg/RData/sim_notsimple/sim2Mle_',k,'_',i,'_',j,'.RData",compress=TRUE)',sep="")
      
      cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",file=paste("sim2Mle_",k,"_",i,"_",j,".R",sep=""))
      ##bashCommands=paste(bashCommands,"nohup R CMD BATCH --vanilla sim_",k,"_",i,"_",j,".R > /dev/null &","\n",sep="")
      ##cat(bashCommands, file="simjobssub.sh")
      if(subjob){
        pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
        pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
        pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
        pbsCommands=paste(pbsCommands,"\n","#$ -N Sim2Mle_",k,"_",i,"_",j,sep="")
        pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH sim2Mle_",k,"_",i,"_",j,".R",sep="")
        cat(pbsCommands,file=paste('runSim2Mle_',k,'_',i,'_',j,'.sh',sep=""),append=FALSE)
        system(paste('chmod u+x runSim2Mle_',k,'_',i,'_',j,'.sh',sep=""))
        system(paste('/opt/sge/bin/lx24-amd64/qsub runSim2Mle_',k,'_',i,'_',j,'.sh',sep=""))
      }
    }
  }
}
