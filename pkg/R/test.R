# rokasPhi <- read.csv("~/proteinevoutk20/pkg/Data/Rokas/genePhi.csv",header=TRUE)
# phi <- read.csv("~/proteinevoutk20/pkg/Data/phi.csv")
# # convert class of data frame from factor to character
# phi <- data.frame(lapply(phi,as.character),stringsAsFactors=FALSE)
# rokasPhi <- data.frame(lapply(rokasPhi,as.character),stringsAsFactors=FALSE)
get_seq <- function(siminfo){
  nsteps <- nrow(siminfo$sim)
  nsites <- ncol(siminfo$sim) - 2
  return(as.numeric(siminfo$sim[nsteps,1:nsites]))
}

gphi <- NULL
pdiff_emp_emp <- NULL
pdiff_new_emp <- NULL
pdiff_emp_new <- NULL
pdiff_new_new <- NULL
dis_diff_emp_emp <- NULL
dis_diff_new_emp <- NULL
dis_diff_emp_new <- NULL
dis_diff_new_new <- NULL
for(gene in 1:5){
  RDatafile <- paste("~/BackupProEvo/Newton/rokas/prunetree/simEqm/gene",gene,".RData",sep="")
  if(file.exists(RDatafile))
    load(RDatafile)
  else{
    cat("gene",gene,", Data file does not exist!")
    gphi <- c(gphi,NA)
    pdiff_emp <- c(pdiff_emp,NA)
    pdiff_new <- c(pdiff_new,NA)
    dis_diff_emp <- c(dis_diff_emp,NA)
    dis_diff_new <- c(dis_diff_new,NA)
  }
  sim_seq <- vector(mode="list",length=4)
  names(sim_seq) <- names(sim)
  for(j in 1:4){
    sim_seq[[j]] <- sapply(sim_info[[j]],get_seq) #matrix with 200 columns
  }
  nsites <- length(index_p)
  gphi <- c(gphi,s)
  pdiff_emp_emp <- c(pdiff_emp_emp,mean(sapply(1:nsim,function(x) sum(sim_seq$emp_emp[,x] != obs.seq)/nsites)))
  pdiff_emp_new <- c(pdiff_emp_new,mean(sapply(1:nsim,function(x) sum(sim_seq$emp_new[,x] != obs.seq)/nsites)))
  pdiff_new_emp <- c(pdiff_new_emp,mean(sapply(1:nsim,function(x) sum(sim_seq$new_emp[,x] != obs.seq)/nsites)))
  pdiff_new_new <- c(pdiff_new_new,mean(sapply(1:nsim,function(x) sum(sim_seq$new_new[,x] != obs.seq)/nsites)))
  
  dis_diff_emp_emp <- c(dis_diff_emp_emp,mean(sapply(1:nsim,function(x) -tail(sim_info$emp_emp[[x]]$dis,1))))
  dis_diff_emp_new <- c(dis_diff_emp_new,mean(sapply(1:nsim,function(x) -tail(sim_info$emp_new[[x]]$dis,1))))
  dis_diff_new_emp <- c(dis_diff_new_emp,mean(sapply(1:nsim,function(x) -tail(sim_info$new_emp[[x]]$dis,1))))
  dis_diff_new_new <- c(dis_diff_new_new,mean(sapply(1:nsim,function(x) -tail(sim_info$new_new[[x]]$dis,1))))
}
save(list=c("gphi",ls()[grep("diff",ls())]),file="adequacy.RData",compress=TRUE)