get.root <- function(mat){
  l <- dim(mat)[2]
  sapply(1:l,function(x) which(mat[,x]==1))
}
simulationQ1 <- function(protein,t,Q=NULL,bf=NULL){
  if(is.null(bf)) bf <- rep(1/20,20)
  if(is.null(Q)) Q <- rep(1,190)
  if(is.matrix(Q)) Qmat <- Q
  else if(is.vector(Q))  
    Qmat <- mat_form_lowtriQ(Q,bf,byrow=TRUE) #transition rate matrix
  l <- length(protein) #number of sites
  P <- expm(Qmat*t)
  res <- protein
  for(i in 1:20){
    index <- which(protein==i)
    res[index] <- sample(1:20,length(index),replace=TRUE,prob=P[i,])
  }
  return(res)
}

find.Q <- function(seq1,seq2,Q=rep(1,6),t=1,method="COBYLA",maxeval="1000",print_level=0){
  fn <- function(Q){
    Q <- exp(Q)
    mumat <- aa_MuMat_form(Q)
    P <- expm(mumat*t)
    prob <- -log(sapply(1:length(seq1), function(i) P[seq1[i],seq2[i]]))
    cat(Q,sum(prob),"\n")
    return(sum(prob))
  }
  lower=rep(-Inf,6)
  upper=rep(Inf,6)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=log(Q),eval_f=fn, lb=lower,ub=upper,opts=opts)
  #optimx(log(Q),fn,hessian=FALSE,method=method,itnmax=itnmax)
  return(res)
}

find.Q1 <- function(seq1,seq2,q=1,t=1,interval=c(-10,10)){
  fn <- function(q){
    q = exp(q)
    Q <- rep(q,6)
    mumat <- aa_MuMat_form(Q)
    P <- expm(mumat*t)
    prob <- -log(sapply(1:length(seq1), function(i) P[seq1[i],seq2[i]]))
    cat(q,sum(prob),"\n")
    return(sum(prob))
  }
  optimize(f=fn,interval=interval)
}

optim.ultrametric <- function(tree,data,model="WAG"){
  
}

internal.br <- function(tree){
  nTips <- length(tree$tip.label)
  edge <- tree$edge
  internal.ind <- which(!(edge[,2] %in% c(1:nTips)))
}