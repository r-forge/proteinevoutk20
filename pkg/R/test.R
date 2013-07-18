num.tip <- 8
setwd(paste("~/BackupProEvo/Newton/balBobyqa_eqm/tip",num.tip,sep=""))
load(paste(num.tip,"tip3000char.RData",sep=""))
opaa <- opaa[index]
dif <- which(opaa!=opaa1)
opaa.dif <- opaa[dif]
opaa1.dif <- opaa1[dif]

dif.data.AA <- matrix(AA[datanum[,dif]],nrow=num.tip)
rownames(dif.data.AA) <- rownames(datanum)
dif.data <- phyDat(dif.data.AA,type="AA")
Qall <- res.op$Qall
tree <- res.op$tree
dif.ind <- attr(dif.data,"index")
opaa_cmp <- function(y)
  sapply(c(opaa[dif[y]],opaa1[dif[y]]),
         function(x) exp(as.vector(ll3m(dif.data,tree,Q=Qall[[x]])$re))[dif.ind[y]])
#################################################################
opaas <- cbind(opaa1[dif],opaa[dif]) #only the cases where inferred opaa is different from true values
#opaas.unique <- cbind(opaa1,opaa)
opaas.unique <- unique.array(opaas) #unique rows
opaas.unique <- opaas.unique[order(opaas.unique[,1]),] #order by the true values
counts <- NULL #counts of each mismatch pair
for(i in 1:dim(opaas.unique)[1]){
  count <- length(intersect(which(opaas[,1]==opaas.unique[i,1]),which(opaas[,2]==opaas.unique[i,2])))
  counts <- c(counts,count)
}
#bubble plot 
symbols(x=opaas.unique[,1],y=opaas.unique[,2],circles=counts,inches=1/4,bg="blue",
        main="mismatched opaas",xlab="start opaa",ylab="inferred opaa")
abline(v=1:20,h=1:20)
#################################################################
opaas.all <- cbind(opaa1,opaa)
opaas.all.unique <- unique.array(opaas.all) #unique rows
opaas.all.unique <- opaas.all.unique[order(opaas.all.unique[,1]),] #order by the true values
all.counts <- NULL #counts of each mismatch pair
for(i in 1:dim(opaas.all.unique)[1]){
  count <- length(intersect(which(opaas.all[,1]==opaas.all.unique[i,1]),which(opaas.all[,2]==opaas.all.unique[i,2])))
  all.counts <- c(all.counts,count)
}
#bubble plot 
symbols(x=opaas.all.unique[,1],y=opaas.all.unique[,2],circles=all.counts,inches=1/4,bg="blue",
        main="Plot of opaas",xlab="start opaa",ylab="inferred opaa")
abline(v=1:20,h=1:20)
#################################################################
##bubble plot:
plot.bubble <- function(opaa,Qall,inches=1/3,grid=TRUE){
  sub1 <- Qall[[opaa]]
  nonzero.ind <- which(sub1>0,arr.ind=T)
  symbols(x=nonzero.ind[,1],y=nonzero.ind[,2],circles=sub1[nonzero.ind],inches=inches,xlim=c(1,20),
          xlab="from", ylab="to",main=paste("opaa = ", opaa))
  if(grid)
    abline(h=1:20,v=1:20)
}
###########################################################################
## look at the nonzero entries in the AArate matrix, and try to figure out the 
## influence of NU_rates on AA_rates
# nonzeroPair <- indAApair[!sapply(indAApair,FUN=is.null)]
# sapply(1:length(nonzeroPair),function(x) nonzeroPair[[x]]$rates)
# sapply(1:length(nonzeroPair),function(x) nonzeroPair[[x]]$tot)
# NU_to_AA <- matrix(0,nrow=6,ncol=length(nonzeroPair))
# for(i in 1:length(nonzeroPair)){
#   rates <- nonzeroPair[[i]]$rates
#   NU_to_AA[,i] <- sapply(1:6,function(x) sum(rates==x))/nonzeroPair[[i]]$tot
# }
# get.root <- function(mat){
#   l <- dim(mat)[2]
#   sapply(1:l,function(x) which(mat[,x]==1))
# }
optimQ <- function(tree,data,Q=rep(1,6),method="SBPLX",maxeval="100",print_level=0, ...){
  Q = Q/Q[6] #make last rate 1
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,Q=Q,...)
  #these don't change with the change of s
  fixmatall = res.initial$fixmatall
  
  ab <- Q[1:5] # optimize on the first 5 rates
  fn = function(ab,tree,data){
    cat(ab,"\n")
    result = -mllm1(data,tree,Q=c(ab,1),fixmatall=fixmatall, ...)$ll$loglik
    return(result)
  }
  lower=rep(0,5)
  upper=rep(Inf,5)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  res$solution = c(res$solution,1) # append the last rate (1) to the rate vector
  return(res)
}

hessianQ <- function(mle,data,tree,...){
  fn <- function(Q){
    print(Q)
    return(-mllm1(data=data,tree=tree,Q=Q,...)$ll$loglik)
  }
  return(hessian(fn,mle,method="Richardson"))
}