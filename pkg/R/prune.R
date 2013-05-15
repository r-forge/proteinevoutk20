source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
######################################################################
# find the best model from prottest result file
get_best_model <- function(filename){
  cm1 <- paste("cat ", filename, " | grep 'Best model according to AIC' |awk '{print $NF}'")
  bmodel <- system(cm1,intern=TRUE)
  bmodel_p <- strsplit(bmodel,"+",fixed=TRUE)[[1]]
  is.G <- "G" %in% bmodel_p #gamma?
  is.F <- "F" %in% bmodel_p #empirical frequencies?
  is.I <- "I" %in% bmodel_p #invariant sites included?
  best.inv = 0
  best.shape = 1
  if(is.G){
    cm <- paste("cat ",filename," | grep -A 4 ","' : ",bmodel,"' | grep shape | awk '{print $NF}'",sep="")
    best.shape = as.numeric(system(cm,intern=TRUE))
  }
  if(is.I){
    cm <- paste("cat ",filename, " | grep -A 4 ","' : ",bmodel,"' | grep invariable | awk '{print $NF}'",sep="")
    best.inv = as.numeric(system(cm,intern=TRUE))
  }
  cm <- paste("cat ",filename, " | grep -A 5 ","' : ",bmodel,"' | grep 'lnL' | awk '{print $NF}'",sep="")
  lnL = as.numeric(system(cm,intern=TRUE))
  ## JTT LG DCMut MtREV MtMam MtArt Dayhoff WAG RtREV CpREV Blosum62 VT HIVb HIVw FLU
  ## "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"
  model <- bmodel_p[1]
  if(model=="MtREV")
    bmodel_p[1] <- "mtREV24"
  else if(model == "MtMam")
    bmodel_p[1] <- "mtmam"
  else if(model == "MtArt")
    bmodel_p[1] <- "mtArt"
  else if(model == "CpREV")
    bmodel_p[1] <- "cpREV"
  else if(model == "MtMam")
    bmodel_p[1] <- "mtmam"
  else if(model == "Blosum62")
    bmodel_p[1] <- "Blossum62"
  return(list(shape=best.shape[1],inv=best.inv[1],lnL=lnL[1],model=bmodel_p))
}

##################################################################
## optimize two branch lengths, one is newly added, the other is one of the splitting 2 branches
optim.br.add <- function(data,tree,new.ext,new.splits,sumsplits,el=NULL,method="COBYLA",maxeval="100", print_level=0, ...){
  if(is.null(el))
  {el <- tree$edge.length}
  ab <- c(el[new.ext],min(sumsplits,el[new.splits[1]]))
  
  res.initial = mllm1(data=data,tree=tree,...)
  #these don't change with the change of s
  Qall = res.initial$Qall
  bfaa = res.initial$bfaa
  #ab <- tree$edge.length[c(new.ext,new.splits[1])] #what if the second number is bigger than sumsplits?
  fn = function(ab,data,tree){
    tree$edge.length[new.ext] = ab[1]
    tree$edge.length[new.splits[1]] = ab[2]
    tree$edge.length[new.splits[2]] = sumsplits - ab[2]
    #print(c(ab,sumsplits-ab[2]))
    #print(tree$edge.length)
    result = -mllm1(data,tree, Qall=Qall,bfaa=bfaa,...)$ll$loglik
    return(result)
  }
  lower=rep(0,2)
  upper=c(Inf,sumsplits)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,
               "stopval"=-Inf,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  #print(res)
  return(res)
}
optim.br.add.emp <- function(data,tree,new.ext,new.splits,sumsplits,el=NULL,method="COBYLA",maxeval="100", print_level=0, ...){
  if(is.null(el))
  {el <- tree$edge.length}
  ab <- c(el[new.ext],min(sumsplits,el[new.splits[1]]))
  fn = function(ab,data,tree){
    tree$edge.length[new.ext] = ab[1]
    tree$edge.length[new.splits[1]] = ab[2]
    tree$edge.length[new.splits[2]] = sumsplits - ab[2]
    #print(tree$edge.length)
    result = -pml(tree=tree,data=data,...)$logLik
    return(result)
  }
  lower=rep(0,2)
  upper=c(Inf,sumsplits)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,
               "stopval"=-Inf,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  #print(res)
  return(res)
}
######################################################################
## prune one branch and grow it back
######################################################################
## analysis under new model
## filename: file that contains the fasta data of nuleotide data
## dtip: the tip to delete from the tree, could be the order or the actual name
## tree: the tree to start with, with all the tips
## ancestral: method of specifying the root states
prune_new <- function(filename,dtip,tree,ancestral="eqm",range=NULL){
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") #make sure it's pruningwise order
    tree <- reorder.phylo(tree,order="pruningwise") #ape function
  data = conv(filename,range=range,type="phyDat")
  if(is.numeric(dtip)) dtip = tree$tip.label[dtip]
  tree_p = drop.tip(tree,dtip) # new tree, after trim the tip
  tree_p <- reorder.phylo(tree_p,order="pruningwise") 
  ##sometimes tree has many branches with same lengths, in that case the comp.tree won't work :(
  tree1 <- tree
  tree1$edge.length <- runif(length(tree$edge.length))
  tree_p1 = drop.tip(tree1,dtip) # new tree, after trim the tip
  tree_p1 <- reorder.phylo(tree_p1,order="pruningwise") 
  
  brs <- comp.tree(tree1,tree_p1)

  br.split <- brs$br.split #branch in tree_p that's splitted into 2
  new.ext <- brs$new.ext #branch in tree that leads to dtip
  new.splits <- brs$new.splits #branches in tree that add to br.split in tree_p
  shareind <- brs$shareind
  
  data_p = subset(data,subset=tree_p$tip.label) # pruned data
  #data_p = phyDat(as.character(data_p),type="AA") #re-weight pruned data
  res_p <- mllm1(data_p,tree_p,s=1,beta=be,gamma=ga,Q=NU_VEC,ancestral=ancestral) #ll of pruned data
  if(test){
    res_op <- res_p
    iter="10"
  }
  else{
    res_op <- optim.mllm1(res_p,optQ=T,optBranch=T,optsWeight=T,optOpw=FALSE,
                          control=list(epsilon=1e-08,hmaxit=50,htrace=1,print_level=0,maxeval="100"),ancestral=ancestral)
    iter="200"
  }
  ## MLE for the parameters
  s = res_op$s
  beta=res_op$GMweights[2]
  gamma = res_op$GMweights[3]
  ancestral = res_op$ll$ancestral
  Q = res_op$Q
  tree_p = res_op$tree
  index_p = attr(data_p,"index")
  opaa_p = res_op$ll$opaa
  bfaa = res_op$bfaa
  dismat = res_op$dismat
  mumat = res_op$mumat
  fixmatall = res_op$fixmatall
  Qall <- res_op$Qall
  
  tree$edge.length[shareind[,"share"]] <- tree_p$edge.length[shareind[,"pshare"]]
  sumsplits = tree_p$edge.length[br.split]
  ## optimize branch lengths on 8-tip tree with parameters just found, and 8-tip data
  br <- optim.br.add(data,tree,new.ext=new.ext,new.splits=new.splits,sumsplits=sumsplits,
                     maxeval=iter,print_level=1,mumat=mumat,fixmatall=fixmatall,ancestral=ancestral,
                     opaa=opaa_p)

  tree$edge.length[new.ext] <- br$solution[1] #assign the tree optimized branch lengths
  tree$edge.length[new.splits[1]] <- br$solution[2]
  tree$edge.length[new.splits[2]] <- sumsplits - br$solution[2]
  
  brlen <- br$solution[1]  # length of the re-grafted branch
  tree1 <- tree #tree with regrafted branch length 0
  tree1$edge.length[new.ext] <- 0 #make that branch 0
  
  nr <- attr(data_p,"nr") #number of different patterns in 7-tip data
  # for site i, find the probability (likelihood) of all 20 states at that site
  # for every site, set the observed state to be i, and then loop through all the 20 states
  state_i <- function(x){
    datai <- data_p #pruned data
    datai$add <- as.integer(rep(x,nr)) #add the pruned tip back
    names(datai)[length(datai)] <- dtip #change the name back
    ll_i <- mllm1(datai,tree1,Qall=Qall,opaa=opaa_p,bfaa=bfaa,ancestral=ancestral)$ll$sitelik #ll for all sites
    return(ll_i)
  }
  sitell <- sapply(1:20,state_i) #loop through all states and put results together in a matrix
  siteprob <- exp(sitell)
  return(list(brlen=brlen,tree=tree1,res=res_op,prob=siteprob))
}

######################################################################
## analysis under empirical model
######################################################################
prune_emp <- function(filename,dtip,tree,model,range=NULL){
  tree <- unroot(tree) #unroot the tree, (unrooted trees are used for empirical models)
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") #make sure it's pruningwise order
    tree <- reorder.phylo(tree,order="pruningwise") #ape function
  data = conv(filename,range=range,type="phyDat")
  if(is.numeric(dtip)) dtip = tree$tip.label[dtip]
  tree_p = unroot(drop.tip(tree,dtip)) # new tree, after trim the tip
  tree_p <- reorder.phylo(tree_p,order="pruningwise") 
  #sometimes tree has many branches with same lengths, in that case the comp.tree won't work :(
  tree1 <- tree
  tree1$edge.length <- runif(length(tree$edge.length))
  tree_p1 = unroot(drop.tip(tree1,dtip)) # new tree, after trim the tip
  tree_p1 <- reorder.phylo(tree_p1,order="pruningwise") 
  
  brs <- comp.tree(tree1,tree_p1)
  br.split <- brs$br.split #branch in tree_p that's splitted into 2
  new.ext <- brs$new.ext #branch in tree that leads to dtip
  new.splits <- brs$new.splits #branches in tree that add to br.split in tree_p
  shareind <- brs$shareind
  
  data_p = subset(data,subset=tree_p$tip.label) # pruned data
  #data_p = phyDat(as.character(data_p),type="AA") #re-weight pruned data
    
  is.G <- "G" %in% model #gamma?
  is.F <- "F" %in% model #empirical frequencies?
  is.I <- "I" %in% model #invariant sites included?
  bf = NULL
  if(is.F) bf = findBf2(data_p)
  k = 1
  shape = 1
  if(is.G) k=4
  res_p <- pml(tree=tree_p,data=data_p,bf=bf,inv=0,k=k,shape=1,model=model[1]) #ll of pruned data
  res_op <- optim.pml(res_p,optGamma=is.G,optInv=is.I,optEdge=TRUE) #optimize parameters
  tree_p <- res_op$tree
  
  tree$edge.length[shareind[,"share"]] <- tree_p$edge.length[shareind[,"pshare"]]
  sumsplits = tree_p$edge.length[br.split]
  #browser()
  br <- optim.br.add.emp(data=data,tree=tree,new.ext=new.ext,new.splits=new.splits,sumsplits=sumsplits,print_level=1,
                         bf=res_op$bf,k=k,shape=res_op$shape,inv=res_op$inv,model=model[1])

  tree$edge.length[new.ext] <- br$solution[1] #assign the tree optimized branch lengths
  tree$edge.length[new.splits[1]] <- br$solution[2]
  tree$edge.length[new.splits[2]] <- sumsplits - br$solution[2]
  
  brlen <- br$solution[1]  # length of the re-grafted branch
  tree1 <- tree #tree with regrafted branch length 0
  tree1$edge.length[new.ext] <- 0 #make that branch 0
  
  nr <- attr(data_p,"nr") #number of different patterns in 7-tip data
  # for site i, find the probability (likelihood) of all 20 states at that site
  # for every site, set the observed state to be i, and then loop through all the 20 states
  state_i <- function(x){
    tip_ind <- which(names(data)==dtip)
    datai <- data #pruned data
    datai[[tip_ind]] <- as.integer(rep(x,nr))
    res <- pml(tree1,datai)
    res <- update(res,bf=res_op$bf,k=k,shape=res_op$shape,inv=res_op$inv,model=model[1]) #ll for all sites
    return(res$siteLik)
  }
  browser()
  sitell <- sapply(1:20,state_i) #loop through all states and put results together in a matrix
  siteprob <- exp(sitell)
  return(list(brlen=brlen,tree=tree1,res=res_op,prob=siteprob))
}
get.brlen <- function(tree,tip){
  br.index <- which(tree$edge[,2]==which(tree$tip.label==tip)) #index of edge.length of the pruned branch
  brlen <- tree$edge.length[br.index]  # length of the re-grafted branch
  return(list(index=br.index,brlen=brlen))
}

comp.tree <- function(tree,ptree){
  dtip <- setdiff(tree$tip.label,ptree$tip.label)
  dtip.ind <- which(tree$tip.label==dtip)
  br.split <- which(!ptree$edge.length %in% tree$edge.length) #branch index in ptree that's splitted into 2
  new.br <- which(!tree$edge.length %in% ptree$edge.length) #newly added branches in bigger tree
  new.br.tos <- tree$edge[new.br,2] #ends of newely added branches
  new.ext <- new.br[new.br.tos==dtip.ind] #newly added branch that leads to the new tip
  new.splits <- new.br[new.br.tos!=dtip.ind] #newly added two branches that add to the old branch
  pshare <- which(ptree$edge.length %in% tree$edge.length) #indices of shared branches in ptree
  share1 <- which(tree$edge.length %in% ptree$edge.length) #indices of shared branches in tree (probably not matched with pshare)
  share <- share1
  for(i in 1:length(pshare)){
    share[i] <- share1[tree$edge.length[share1]==ptree$edge.length[pshare[i]]]
  }
  shareind <- cbind(pshare,share)
  if(sum(tree$edge.length[new.splits])!=ptree$edge.length[br.split])
    stop("something is wrong!")
  return(list(br.split=br.split,new.ext=new.ext,new.splits=new.splits,shareind=shareind))
}
######################################################################
## this version's precondition: bigger tree and the tip that's deleted
## instead of a  bigger tree and a smaller tree
comp.tree1 <- function(tree,dtip){
  rooted <- is.rooted(tree)
  tree$edge.length <- runif(length(tree$edge.length))
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") #make sure it's pruningwise order
    tree <- reorder.phylo(tree,order="pruningwise") #ape function
  if(is.numeric(dtip)) dtip = tree$tip.label[dtip]
  ptree = drop.tip(tree,dtip) # new tree, after trim the tip
  ptree <- reorder.phylo(ptree,order="pruningwise")
  if(!rooted) ptree <- unroot(ptree)
  #sometimes tree has many branches with same lengths, in that case the comp.tree won't work :(
  dtip.ind <- which(tree$tip.label==dtip)
  br.split <- which(!ptree$edge.length %in% tree$edge.length) #branch index in ptree that's splitted into 2
  new.br <- which(!tree$edge.length %in% ptree$edge.length) #newly added branches in bigger tree
  new.br.tos <- tree$edge[new.br,2] #ends of newely added branches
  new.ext <- new.br[new.br.tos==dtip.ind] #newly added branch that leads to the new tip
  new.splits <- new.br[new.br.tos!=dtip.ind] #newly added two branches that add to the old branch
  pshare <- which(ptree$edge.length %in% tree$edge.length) #indices of shared branches in ptree
  share1 <- which(tree$edge.length %in% ptree$edge.length) #indices of shared branches in tree (probably not matched with pshare)
  share <- share1
  for(i in 1:length(pshare)){
    share[i] <- share1[tree$edge.length[share1]==ptree$edge.length[pshare[i]]]
  }
  shareind <- cbind(pshare,share)
  if(sum(tree$edge.length[new.splits])!=ptree$edge.length[br.split])
    stop("something is wrong!")
  return(list(br.split=br.split,new.ext=new.ext,new.splits=new.splits,shareind=shareind))
}
