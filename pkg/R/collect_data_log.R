##numer of sites in the simulation
#numsites <- c(1000,2000,3000)
numsites <- c(4000,6000,10000,20000)
##s values used in simulation
svals <- c(0.1,0.5,0.9)
##tips of trees used in simulation
#tip <- c(2,4,8,16)
## edge length scale factor
#esvalues <- c(1,10)
##result to store all the data
#opaa <- c(1:5)
datares <- array(list(NULL),dim=c(3,3,5))
##change this to "simMle" for information on MLE without log
fileprefix <- "sim_"
DataFile <- "summary.RData"


  for(sitesindex in 1:3){
    for(sindex in 1:3){
      for(opaa in 1:5){
        cat(numsites[sitesindex], "sites, s=",svals[sindex],"optimal amino acid: ",opaa, "\n")
        RDataName = paste(fileprefix,numsites[sitesindex],"_",svals[sindex],"_",opaa,".RData",sep="")
        if(file.exists(RDataName)){
          load(RDataName)
          datares[[sitesindex,sindex,opaa]]$op <- sw_op$par
          datares[[sitesindex,sindex,opaa]]$est <- sw_est$par
          datares[[sitesindex,sindex,opaa]]$opw <- c(sw_opw$s,sw_opw$GMweights[2:3])
          datares[[sitesindex,sindex,opaa]]$weights <- sw_opw$opw
         }
       else
         cat("file doesn't exist!","\n")
       }
    }
  }

save(datares,file=DataFile)

getMle <- function(sind,opaa,paraind){
  res <- matrix(0,4,3)
  res[1,] <- rep(svals[sind],3)
  for(i in 1:3){
    res[2:4,i] = c(datares[[i,sind,opaa]]$est[paraind],datares[[i,sind,opaa]]$op[paraind],datares[[i,sind,opaa]]$opw[paraind])
  }
  return(res)
}
getMle.s <- function(opaa){
  cbind(getMle(1,opaa,1),getMle(2,opaa,1),getMle(3,opaa,1))
}
getMle.beta <- function(opaa){
  cbind(getMle(1,opaa,2),getMle(2,opaa,2),getMle(3,opaa,2))
}
getMle.gamma <- function(opaa){
  cbind(getMle(1,opaa,3),getMle(2,opaa,3),getMle(3,opaa,3))
}
## ML estimates of s and weights, for all genes in Rokas's data
rokasres = matrix(0,106,8)
Qres = matrix(0,106,6)
elres = matrix(0,106,14)
for(geneindex in 1:106){
  cat("gene", geneindex, "\n")
  RDataName = paste("gene",geneindex,"_s_weight.RData",sep="")
  if(file.exists(RDataName)){
    load(RDataName)
    rokasres[geneindex,] = c(res_sweight$par,res_sweight$value,res_op$s,res_op$GMweights[2:3],res_op$ll$loglik)
    Qres[geneindex,] = res_op$Q
    elres[geneindex,] = res_op$tree$edge.length
  }
  else
    cat("file doesn't exist!","\n")
}