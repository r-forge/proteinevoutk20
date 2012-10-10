##numer of sites in the simulation
numsites <- c(500,800,1000,2000,3000)
##s values used in simulation
svals <- c(0.1,0.3,0.5,0.7,0.9)
##tips of trees used in simulation
tip <- c(2,4,8,16)
## edge length scale factor
esvalues <- c(1,10)
##result to store all the data
datares <- array(list(NULL),dim=c(4,5,5,2))
##change this to "simMle" for information on MLE without log
fileprefix <- "sim_"
DataFile <- "summary.RData"

for(treeindex in 1:4){
  for(sitesindex in 1:5){
    for(sindex in 1:5){
      for(esind in 1:2){
        cat("tree", treeindex, ",", numsites[sitesindex], "sites, s=",svals[sindex],"edge length *",esvalues[esind], "\n")
        RDataName = paste(fileprefix,treeindex,"_",numsites[sitesindex],"_",svals[sindex],"_",esvalues[esind],".RData",sep="")
        if(file.exists(RDataName)){
          load(RDataName)
          datares[[treeindex,sitesindex,sindex,esind]]$op <- sw_op$par
          datares[[treeindex,sitesindex,sindex,esind]]$est <- sw_est$par
          datares[[treeindex,sitesindex,sindex,esind]]$mode <- sw_mode$par
         }
       else
         cat("file doesn't exist!","\n")
       }
    }
  }
}
save(datares,file=DataFile)

