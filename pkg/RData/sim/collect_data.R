##numer of sites in the simulation
numsites <- c(100,300,500,700,1000)
##s values used in simulation
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
##tips of trees used in simulation
tip <- c(4,4,8,8,12,12)
##result to store all the data
datares <- array(list(NULL),dim=c(6,5,6))
##change this to "simMle" for information on MLE without log
##fileprefix <- "simMleLog"
##DataFile <- "collect_data_log.RData"
fileprefix <- "simMle"
DataFile <- "collect_data.RData"
for(treeindex in 1:6){
  for(sitesindex in 1:5){
    for(sindex in 1:6){
      print(paste("tree ",treeindex,", ",numsites[sitesindex], " sites, s = ",svalues[sindex], sep=""))
      RDataName = paste(fileprefix,treeindex,"_",numsites[sitesindex],"_",svalues[sindex],".RData",sep="")
      if(file.exists(RDataName)){
        print("file exists!")
        load(RDataName)
        datares[[treeindex,sitesindex,sindex]]$para <- log(c(svalues[sindex],be,ga))
        datares[[treeindex,sitesindex,sindex]]$mle <- mle$par
        rm(mle)
        if(exists("hes")){
          datares[[treeindex,sitesindex,sindex]]$hessian <- hes
          datares[[treeindex,sitesindex,sindex]]$eigenval <- eigen(hes)$values
          datares[[treeindex,sitesindex,sindex]]$invhes <- solve(hes)
          print("remove hessian now")
          rm(hes)
        }
        else{
          print("hessian not calculated")
##           datares[[treeindex,sitesindex,sindex]]$hessian <- NULL
##           datares[[treeindex,sitesindex,sindex]]$eigenval <- NULL
##           datares[[treeindex,sitesindex,sindex]]$invhes <- NULL
        }
      }
      else
        print("file doesn't exist!")
    }
  }
}
save(datares,file=DataFile)

