##################################################################
##   load the required libraries ##
##################################################################
library(seqinr) #package to convert codons to amino acids
library(phangorn) #phylogenetics package
library(optimx) #optimization package
library(expm) #matrix exponentiation
#library(parallel) #parallel computing
library(multicore)
library(minqa) #optimization function bobyqa
library(mgcv) #find the unique rows in a matrix - uniquecombs
library(numDeriv) # calculate hessian of a function at a point
##################################################################
#directory of the data need to read, including gene, tree, Grantham data
datadir <- "~/proteinevoutk20/pkg/Data/"
##################################################################
##   Read in parameters and data  ##
##################################################################
##Nucleotide mutation rates estimated from Rokas's data (using PAUP):
##Model GTR, with Gamma rate distribution and invariable sites and etc
##Nuvec <-  c(2.94194,8.23705,1.39133,2.33991,14.86435,1.00000)
##Model GTR
Nu_vec <- c(1.96575,4.08655,1.39431,1.46172,6.36024,1.00000)
Nu_vec_JC <- rep(1.0,6)
##empirical base frequencies
freq <- c(0.31065,0.18365,0.20955,0.29615)
##Best tree with branch lengths output by PAUP
tree <- read.nexus(paste(datadir,"GTR.tre",sep=""))
##Read the properties of amino acids, including c(composition),p(polarity) and 
##v(molecular volume), save the data.frame
GMcpv <- read.csv(paste(datadir,"Grantham_cpv.csv",sep=""),header=TRUE, sep=",",row.names=1)
GMcpv <- data.matrix(GMcpv)

##################################################################
##  Grantham Matrix     ##
##################################################################
##Build the Grantham matrix with coefficients for all components as parameters
##alpha for composition,beta for polarity,gamma for molecular volume
##The mean of the distance is normalized to 1,(mean of only the upper triangle elements of the matrix)
##The parameter values in Grantham matrix are:
##alpha=1.833, beta=0.1018,gamma=0.000399
##First find the avg chemical distance, and then the weight is (1/avg_dis)^2

##this function computes the weights that are used by Grantham
## (dispite the name, they are not the average distances though, instead they are functions of average distances)
avg_dis <- function(vector){
  l <- length(vector)
  Sum <- 0
  ## find the sum of all pairwise chemical distances
  for(i in 1:(l-1)){
    Sum <- Sum + sum(abs(vector[i]-vector[(i+1):l]))
  }
  ct <- l*(l-1)/2 #total combination count
  return((ct/Sum)^2)
}
##Grantham weights
al <- avg_dis(GMcpv[,1])
be <- avg_dis(GMcpv[,2])
ga <- avg_dis(GMcpv[,3])


###Given the data (20 by 3 matrix/data frame) from Grantham on all three factors, and the weights,
###Find the distance matrix for all the amino acids (20 by 20) with mean 1
###Grantham matrix has mean 100 instead of 1
GM_cpv <- function(datamatrix, alpha=al, beta=be, gamma=ga){
  A <- datamatrix
  DistanceMatrix <- matrix(0,nrow=20,ncol=20) #diagonal entries have values 0
  for(i in 1:19){
  for(j in (i+1):20)
  	DistanceMatrix[i,j] <- (alpha*(A[i,1] - A[j,1])^2 + beta*(A[i,2]-A[j,2])^2+gamma*(A[i,3]-A[j,3])^2)^(1/2)
  }
  ##make a symmetric matrix
  result <- DistanceMatrix + t(DistanceMatrix)
  ##normalize the mean to 1
  m <- sum(result)/(400-20) #exclude the diagonal entries
  return(result/m)
}
##Default Grantham matrix (one used in Grantham's paper)
GM <- GM_cpv(GMcpv,al,be,ga)


##################################################################
##   construct matrices for mutation and substitution  ##
##################################################################
##Short names for all the amino acids, in alphabetical order
aa <- s2c("ARNDCQEGHILKMFPSTWYV") #amino aicds in alphabetical order
Nu <- s2c("acgt") #nucleotide ACGT - now the order is changed !!!
stop_cd <- c("taa","tag","tga") #stop codons
nu_list <- NULL #list of all 64 codons (including stop codons)
for(i in 1:4){
  for(j in 1:4){
    for(k in 1:4){
      nu_list <- c(nu_list,c2s(c(Nu[i],Nu[j],Nu[k])))
    }
  }
}
cds <- nu_list[!nu_list %in% stop_cd] #61 non-stop codons

##################################################################
##   Read gene data ##
##################################################################
## read in fasta file of nucleotide, convert it into amino acid data
conv <- function(filename){
  levels <- aa #amino acids in alphabeticla order (not the single letter names, the 3-letter names)
  data <- read.fasta(file=filename) #read fasta file including nucleotide data
  dna.to.aa <- function(x){ #given the order of the data, convert nucleotide data to amino acid data
     aas <- translate(seq=data[[x]])
     aan <- as.numeric(factor(aas,levels)) #convert from letters to numbers, as levels
  }
  seqdata <- t(sapply(1:length(data),dna.to.aa,simplify="array")) #sapply and simplify to array
  dimnames(seqdata)[[1]]<- attr(data,"name") #Change the dimnames to the species names
  return(seqdata)
}

##   Read in 106 gene data ##
l <- 106
data <- vector("list",length=l)
##data stores the data of all genes in rokas's data
system.time(
  for(i in 1:l){
    data[[i]] <- conv(paste(datadir,"gene",i,".fasta",sep=""))
  }
)

## Find the most frequent element of a vector
## It works for numbers and characters
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#############################################################################
##Given a vector of length 6 (number of parameters for GTR matrix)
##generate the GTR rate matrix Q
mat_form <- function(vec){
  rate_mat <- matrix(0,nrow=4,ncol=4)
  rate_mat[1,2] <- vec[1]
  rate_mat[1,3] <- vec[2]
  rate_mat[1,4] <- vec[3]
  rate_mat[2,3] <- vec[4]
  rate_mat[2,4] <- vec[5]
  rate_mat[3,4] <- vec[6]
  rate_mat[lower.tri(rate_mat)] <- rate_mat[upper.tri(rate_mat)] #matrix is symmetric
  diag(rate_mat) <- sapply(1:4,function(x) -sum(rate_mat[x,-x])) #row sums equal to 0
  dimnames(rate_mat)[[1]] <- Nu
  dimnames(rate_mat)[[2]] <- Nu
  return(rate_mat)
}

##Find the scale for the subsitution rate matrix Q
scale.rate <- function(Q){
  Pi <- as.vector(expm(Q*10000)[1,]) #equilibrium frequency vector
  return(-sum(Pi*diag(Q)))
}
#############################################################################
##A function I need to do the matrix exponential correctly
##Check if all the entries in x are finite, if yes, return TRUE, otherwise return FALSE
is.nan.inf <- function(x){
  return(all(is.finite(x)==TRUE))
}
## A modified version of expm that make sure there is no Inf and NaN in the result of matrix exponential
expm.m <- function(x){
   edefault <- expm(x) #expm from expm package with default method "higham08.b"
  if(is.nan.inf(edefault)==TRUE) #if there is no Inf or NaN values in the result, return the result
    return(edefault)
  else{
    #print(x)
    eward77 <- expm(x,method="Ward77") #method Ward77
    if(is.nan.inf(eward77)==TRUE)
      return(eward77)
    else{
      epade <- expm(x,method="R_Pade") #method Pade
      if(is.nan.inf(epade)==TRUE)
        return(epade)
      else{
        etaylor <- expm(x,method="taylor") #method Taylor
        if(is.nan.inf(etaylor)==TRUE)
          return(etaylor)
        else{
          return(eward77)
          warning("all methods give Inf or NaN",immediate.=TRUE)
        }
      }
    }
  }
}

#############################################################################
#Find the mutation rate matrix (20 by 20) from the 6 rates in GTR mutation rate
#matrix for nucleotides
#the rate G<->C is normalized to 1, now the scale is to scale the rates to real values
aa_MuMat_form <- function(vec,scale=T){
  Q <- mat_form(vec) #rate matrix for nucleotides
  if(scale==T){ #scale the rate matrix
    Scale = scale.rate(Q)
    Q=Q/Scale
  }
  ##mutation matrix for 61 codons
  codon_array <- array(0,dim=c(61,61),dimnames=list(cds,cds))
  for(m in 1:61){                       #loop through all 61 codons
    cd_str <- cds[m] #mth codon in string format
    cd <- s2c(cd_str)   #one codon, convert from string to char vector
    ##print(paste("starting from",cd_str,":" ))  #for testing purpose
    for(i in 1:3){           #every nucleotide in the codon can mutate
      ngb <- cd #neighboring codon, set it equal to the starting codon for now
         for(j in 1:4){
        if(Nu[j]!=cd[i]){ #nucleotide has to change to a different one
          ngb[i] <- Nu[j] 
          a_new <- translate(ngb)       #new amino acid
          ngb_str <- c2s(ngb)
          ##print(paste("neighbor of",cd_str,":",ngb_str)) #testing
          if(a_new != "*")
            codon_array[cd_str,ngb_str] <- codon_array[cd_str,ngb_str] + Q[cd[i],Nu[j]]
        }
      }
    }
  }
  ##diagonal equals negative of sum of off-diagonal entries
  diag(codon_array) <- sapply(1:61,function(x) -sum(codon_array[x,-x])) #row sums equal to 0
  ##equilibrium frequencies for 61 codons
  cd_eqlm <- expm.m(codon_array*10000)[1,]
  
  ##20*20 array: mutations rates between amino acids
  arr <- array(0,dim=c(20,20),dimnames=list(aa,aa)) #mutation rate matrix for amino acids
  for(i in 1:61){
    for(j in 1:61){
      from_codon <- cds[i]
      to_codon <- cds[j]
      from_aa <- translate(s2c(from_codon))
      to_aa <- translate(s2c(to_codon))
      arr[from_aa,to_aa] <- arr[from_aa,to_aa] + cd_eqlm[from_codon]*codon_array[from_codon,to_codon]
    }
  }
  dimnames(arr) <- NULL
  return(arr)
}
## Amino acid mutation rate matrix for the GTR model of nucleotide
mumat <- aa_MuMat_form(Nu_vec)
mumat_JC <- aa_MuMat_form(Nu_vec_JC)
#############################################################################
##find the physiochemical distance vector between two proteins, given the distance matrix
pchem_d <- function(protein1, protein2,DisMat){
  if(length(protein1)!=length(protein2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
  site_d <- function(k){
    return(DisMat[protein1[k],protein2[k]])
  }
  d <- sapply(c(1:length(protein1)),site_d,simplify=TRUE) 
  return(d)
}
#################################################################
#Functionality function, arguments: 
#d: vector of distances between the optimal and observed amino acids
#s: seleciton coeffecient, can be passed as a vector or a number if it's the same for all sites
Ftny <- function(d, s){  
  if((length(s)==1)&&(length(d)!=1)){ #if s is given as a scalar, then treat it to be the same across all sites
    s <- rep(s,length(d))
  }
  result <- length(d)/(sum(1+d*s)) #harmonic mean of ftny at all sites (F_i = 1 + d_i * s_i)
  return(result)
} #d and s are given

#Functionality of a protein given the optimal protein (aa sequence) and parameters
Ftny_protein <- function(protein,protein_op,s,DisMat){
  d <- pchem_d(protein,protein_op,DisMat)
  return(Ftny(d,s))
}

#fixation probability, since we consider only site by site, so n=1
#d1, d2: distances (vectors)
#s: selection (vector or scalar)
#Ne: effective population size
#q, Phi: constants
#C: cost function, linear function of the length of the protein
#C_n = a1 + a2*n, cost --  linear function of the length
fix <- function(d1,d2,s,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  if((d1==d2)||(s==0)) #When the fitnesses are the same, neutral case, pure drift
    return(1/(2*Ne))
  else{
    fit_ratio <- exp(-C*Phi*q*s*(d1-d2))
    if(fit_ratio==Inf)
      return(0)
    else if(fit_ratio==1)
      return(1/(2*Ne))
    else
      return((1-fit_ratio)/(1-fit_ratio^(2*Ne)))
  }
}


#fixation probability of a mutant with initial freq 1/(2Ne) in a diploid population
#with proteins as parameteres instead of distance vectors for the proteins
fix_protein <- function(protein1, protein2, protein_op, s, DisMat,
                        C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  d1 <- pchem_d(protein1,protein_op,DisMat)
  d2 <- pchem_d(protein2,protein_op,DisMat)
  return(fix(d1,d2,s,C,Phi,q,Ne))
}

######################################################

#Function to generate the transition matrix (m*m for now) 
#given the optimal amino acid. This is only for one site (amino acid)
#parameter: aa_op: optimal amino acid at this site
#s: selection coefficient(SCALAR, since it's only one site)
#m: number of states (amino acids) considered
#DisMat: distance matrix between amino acids
#MuMat: mutation rate matrix between amino acids
mat_gen_indep <- function(aa_op,s,DisMat,MuMat,C=2, Phi=0.5,q=4e-7,Ne=1.36e7){
  m = 20
  mat <- matrix(0,nrow=m,ncol=m) #set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      mat[i,j] <- 2*Ne*MuMat[i,j]*fix_protein(i,j,aa_op,s,DisMat,C,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- 2*Ne*MuMat[j,i]*fix_protein(j,i,aa_op,s,DisMat,C,Phi,q,Ne) #symmetric entry
    }#end for j
    mat[i,i] <- -(sum(mat[i,])) #diagonal entry; sum of each row is 0
  }#end for i
  mat[m,m] <- -sum(mat[m,])
  return(mat)
}

######################################################
##likelihood function for 1 site
ll_site <- function(tree,data,optimal,s,MuMat,alpha=al, beta=be, gamma=ga,
                    bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
    ##If the given tree is not rooted and binary, then throw error and exit
    if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
    m = 20
    ##if the base frequencies are not specified, do a uniform distribution
    if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
    GM1 = GM_cpv(GMcpv,alpha,beta,gamma)
    Q = mat_gen_indep(optimal,s,GM1,MuMat,C,Phi,q,Ne) #transition rate matrix for the site, given the optimal aa
    Q = scaleQ(Q,bf)
    tree <- ape:::reorder.phylo(tree,"p") #reorder the tree in pruningwise order
    edge = tree$edge #edges
    nNodes = max(edge) #number of nodes in the tree
    probvec = matrix(NA,nNodes,m) #probability of gettting different states at nodes that evolve to the current sequences
    parent <- as.integer(edge[, 1]) #parents of the edges
    child <- as.integer(edge[, 2]) #children of the edges
    root <- as.integer(parent[!match(parent, child, 0)][1])  
    tip <- as.integer(child[!match(child,parent,0)])
    init.tip <- function(x){ #initiate the vector for the tips, 1 for the tip state, 0 otherwise
      vec <- rep(0,m)
      vec[data[x]] <- 1
      vec
    }
    probvec[1:length(tip),] <- t(sapply(1:length(tip),init.tip)) #all tips
    tl = tree$edge.length #lengths of the edges
    for(i in 1:tree$Nnode){ #for each interior node calculate the probability vector of observing 1 of 20 states
        from = parent[2*i] #parents
        to = child[(2*i-1):(2*i)] #direct descendents
        t_left = tl[2*i-1] #left branch length
        t_right = tl[2*i] #right branch length
        v.left <- expm.m(Q*t_left) #probabilities of transition from one state to another after time t
        v.right <- expm.m(Q*t_right)
        probvec[from,] <- as.vector((v.left%*%probvec[to[1],])*(v.right%*%probvec[to[2],])) #pruning, vector form
        check.sum <- sum(probvec[from,])
        if(check.sum==0) #probability is very very low
          warning("numerical overflow",immediate.=TRUE)
    }
    return(as.numeric(probvec[root,] %*% bf))
      #return(list(ll=max(probvec[root,]),root=which.max(probvec[root,]))) #with the corresponding root returned 
      #return(max(probvec[root,])) #just the value
}
#vectorized version of  ll_site, on optimal aa and root aa
ll_site_vec <- Vectorize(ll_site,vectorize.args=c("optimal"))
#############################################################################
## Remove columns with NA in the data
PruneMissing <- function(x){
  ##Find the array indices of the NA entries
  naInd <- which(is.na(x),arr.ind=T)
  dimnames(naInd) <- NULL
  ## column indices
  naCol <- naInd[,2]
  x[,-naCol]
}
#############################################################################
##Likelihood of data on a tree, given the selection coefficient "s" and the optimal amino acid sequence "protein_op"
##If protein_op (optimal protein) is not given, get it from the most frequent amino acids appeared in the data
ll_indep <- function(s,alpha,beta,gamma,MuMat,tree,data,protein_op=NULL,bf=NULL,
                             C=2,Phi=0.5,q=4e-7,Ne=1.36e7,majority_op=T){
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  ##Find the unique columns (site data) and occurences of each unique column
  m = 20
  data.u <- t(uniquecombs(t(data))) # if finds uniqe rows, so use the function on transpose
  ind <- attr(data.u,"index") #corresponding row numbers in unique matrix from original matrix
  occu <- as.numeric(table(ind)) # occurences of each unique column
  if(majority_op) #Use the most frequent amino acid, what if there are more than 1 with highest frequency??
     protein_op <- apply(data.u,2,Mode)
  ll.one <- function(x){ #for each unique data, find ll_site, and multiply the -loglikelihood by the occurance
    if(!is.null(protein_op))
      return(-log(ll_site(tree,data.u[,x],protein_op[x],s,MuMat,alpha,beta,gamma,m,root,bf,C,Phi,q,Ne))*occu[x])
    else
      return(-log(max(ll_site_vec(tree,data.u[,x],1:20,s,MuMat,alpha,beta,gamma,m,root,bf,C,Phi,q,Ne)))*occu[x])
  }
  lls <- unlist(lapply(1:length(occu),ll.one))
  if(all(is.finite(lls)==T))
    sum(lls)
  else
    return(10^7)
  #sum(unlist(mclapply(1:length(occu),ll.one))) #parallel version, DO NOT parallelize at this level!!!
}
#system.time(res3 <- ll_indep(0.1,al,be,ga,mumat,tree,data[[3]]))
#############################################################################
##This function finds the MLE estimator for "s" only, given that all other parameters are known,
##including weights(alpha, beta, gamma), mutation rates, C, q, Phi, and Ne
MLE_GTR <- function(start_pt,lowerb,upperb,tree,data,alpha,beta,gamma,MuMat,m=20,protein_op=NULL,root=NULL,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7,trace=0){
  negloglike <- function(s){ #The function to minimize
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,C,Phi,q,Ne))
  }
  ## Delete columns with NA in data
  if(length(which(is.na(data)))!=0) {
    warning("NA in data, pruning performed")
    data <- PruneMissing(data)
  }
  
  #method nlminb using PORT routines
    ans <- nlminb(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(trace=trace))
    #ans <- nmkb(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(trace=TRUE))
  return(ans)
}
##find the MLE of s on the log scale
MLE_GTR_log<- function(start_pt,tree,data,alpha,beta,gamma,MuMat,m=20,protein_op=NULL,root=NULL,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7,trace=0){
  negloglike <- function(s){ #The function to minimize
    return(ll_indep(exp(s),alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,C,Phi,q,Ne))
  }
  ## Delete columns with NA in data
  if(length(which(is.na(data)))!=0) {
    warning("NA in data, pruning performed")
    data <- PruneMissing(data)
  }
  
  ans <- nlminb(start_pt,negloglike,control=list(trace=trace))
  return(ans)
}
#############################################################################
##Given values for beta and gamma, find the MLE's of s for all 106 genes
MLE.s <- function(x,generange,multicore=FALSE){
  Beta <- x[1]
  Gamma <-  x[2]
  mle.s.one <- function(k){
    print(paste("start optimization on gene ", k, sep=""))
    mle <- MLE_GTR(1,0,1e5,tree,data[[k]],al,Beta,Gamma,mumat)
    print(paste("finish optimization on gene ", k, sep=""))
    return(mle)
  }
  if(multicore)
    mclapply(generange,mle.s.one)
  else
    lapply(generange,mle.s.one)
}

MLE.s_log<- function(x,generange,multicore=FALSE){
  Beta <- exp(x[1])
  Gamma <- exp(x[2])
  mle.s.one <- function(k){
    print(paste("start optimization on gene ", k, sep=""))
    mle <- MLE_GTR_log(-1,tree,data[[k]],al,Beta,Gamma,mumat)
    print(paste("finish optimization on gene ", k, sep=""))
    return(mle)
  }
  if(multicore)
    mclapply(generange,mle.s.one)
  else
    lapply(generange,mle.s.one)
}

#system.time(res <- MLE_GTR(1,0,1e4,tree,data[[2]],al,be,ga,mumat))

### optimization of beta and gamma
### starting point, lower bounds, upper bounds, trace, generange, multicore
MLE.bg <- function(start_pt,lowerb,upperb,generange, trace=0,multicore=FALSE){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  mle.bg.s <- function(x){
    mle <- MLE.s(x,generange,multicore=multicore) #call the previous function to optimize s for all genes
    mle.val <- sapply(1:length(generange),function(ind) mle[[ind]]$objective) # best -loglikelihood values
    return(sum(mle.val)) #summation of all values
  }
  ans <- nlminb(start_pt,mle.bg.s,lower=lowerb,upper=upperb,control=list(trace=trace))
  return(ans)
}

MLE.bg_log <- function(start_pt,generange,trace=0,multicore=FALSE){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  mle.bg.s <- function(x){
    mle <- MLE.s_log(x,generange,multicore=multicore) #call the previous function to optimize s for all genes
    mle.val <- sapply(1:length(generange),function(ind) mle[[ind]]$objective) # best -loglikelihood values
    return(sum(mle.val)) #summation of all values
  }
  ans <- nlminb(start_pt,mle.bg.s,control=list(trace=trace))
  return(ans)
}
##This function finds the MLE estimators for s, and the cpv weights. Here alpha is given, beta and gamma are being
##estimated, since only 2 of them are free parameters.
##Brian suggested that the restriction should be that alpha + beta + gamma = 1, DO THIS LATER!
MLE_sw<- function(start_pt,lowerb,upperb,tree,data,m=20,alpha,MuMat,
                      protein_op=NULL,root=NULL,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.37e07,trace=0){
  negloglike <- function(para){
    s = para[1]
    beta = para[2]
    gamma = para[3]
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,C,Phi,q,Ne))
  }
  #ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,method="nlminb",hessian=TRUE,control=list(trace=trace))
  #ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,hessian=TRUE,control=list(trace=trace,all.methods=TRUE))
  ans <- nlminb(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(trace=trace))
  #ans <- bobyqa(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(iprint=3))
  return(ans)
}
## optimization on the log of parameters, so that constraints won't be needed and no worries about
## negative parameter values
MLE_sw_log<- function(start_pt,tree,data,m=20,alpha,MuMat,
                  protein_op=NULL,root=NULL,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.37e07,trace=0){
  negloglike <- function(para){
    s = exp(para[1])
    beta = exp(para[2])
    gamma = exp(para[3])
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,C,Phi,q,Ne))
  }
  #ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,method="nlminb",hessian=TRUE,control=list(trace=trace))
  #ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,hessian=TRUE,control=list(trace=trace,all.methods=TRUE))
  ans <- nlminb(start_pt,negloglike,control=list(trace=trace))
  #ans <- bobyqa(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(iprint=3))
  return(ans)
}
###############################################################################
### objects used in grid search for the best estimates of beta and gamma
#  l <- 15
#  beta <- seq(-15,0,length.out=(l+1))[-1]
#  gamma <- seq(-15,0,length.out=(l+1))[-1]

# grid.1.beta <- seq(0.6,0.75,length.out=(l+1))[-1]
# grid.1.gamma <- seq(0,0.1,length.out=(l+1))[-1]
# 
# grid.3.beta <- seq(0.9,1.0,length.out=(l+1))[-1]
# grid.3.gamma <- grid.1.gamma
# 
# grid.5.beta <- seq(0.45,0.55,length.out=(l+1))[-1]
# grid.5.gamma <- seq(0.05,0.15,length.out=(l+1))[-1]
# 
# grid.6.beta <- seq(0.80,0.90,length.out=(l+1))[-1]
# grid.6.gamma <- grid.5.gamma
# AAMat1 <- function(CdMat,CdBf){
#   CdBf= CdBf/sum(CdBf)
#   Q = matrix(0,20,20)
#   for(i in 1:19){
#     fcodons = cdlist[[i]]
#     fcfq = CdBf[fcodons] #frequencies of codons coding for this amino acid
#     fcfq = fcfq/sum(fcfq) # normalize the sum to 1
#     for(j in (i+1):20){
#       if(i!=j){
#         tcodons = cdlist[[j]]
#         tcfq = CdBf[tcodons] #frequencies of codons coding for this amino acid
#         tcfq = tcfq/sum(tcfq) # normalize the sum to 1
#         Q[i,j] = sum(diag(fcfq) %*% CdMat[fcodons,tcodons])
#         Q[j,i] = sum(diag(tcfq) %*% CdMat[tcodons,fcodons])
#       }
#     }
#   }
#   diag(Q) = -rowSums(Q)
#   return(Q)
# }
# aa_MuMat_form1 <- function(vec=rep(1,6),bf=rep(0.25,4)){
#   Q <- mat_form_lowtriQ(vec,bf) #rate matrix for nucleotides    
#   dimnames(Q) = list(Nu,Nu)
#   ##mutation matrix for 61 codons
#   codon_array <- array(0,dim=c(61,61),dimnames=list(CDS,CDS))
#   for(m in 1:61){                       #loop through all 61 codons
#     cd_str <- CDS[m] #mth codon in string format
#     cd <- s2c(cd_str)   #one codon, convert from string to char vector
#     for(i in 1:3){           #every nucleotide in the codon can mutate
#       ngb <- cd #neighboring codon, set it equal to the starting codon for now
#       for(j in 1:4){
#         if(Nu[j]!=cd[i]){ #nucleotide has to change to a different one
#           ngb[i] <- Nu[j] # change the ith position in codon to jth nucleotide
#           a_new <- translate(ngb)       #new amino acid
#           ngb_str <- c2s(ngb)
#           if(a_new != "*")
#             codon_array[cd_str,ngb_str] <- codon_array[cd_str,ngb_str] + Q[cd[i],Nu[j]]
#         }
#       }
#     }
#   }
#   diag(codon_array) <- -rowSums(codon_array) #row sums equal to 0
#   cdbf <- freq_codon(bf)
#   arr <- AAMat(codon_array,cdbf)
#   return(arr)
# }
# given optimal amino acid, find the rate matrix Q and the eigen-decomposition of Q
eigOp <- function(op_aa,s,DisMat,MuMat,bf=rep(1/20,20),C=2,Phi=0.5,q=4e-07,Ne=1.36e07){
  mat = mat_gen_indep(op_aa,s,DisMat,MuMat,C,Phi,q,Ne)
  res = eigQ(t(mat),bf) #in order to use functions (C and R) from phangorn
  return(res)
}
## given s, distance matrix, and mutation rate matrix, find rate matrices for all 20 amino acids as optimal
## And then find their eigen decomposition 
eigAllaa <- function(s,DisMat,MuMat,bf=rep(1/20,20),C=2, Phi=0.5,q=4e-7,Ne=1.36e7){
  res = vector(mode="list",length=20)
  for(i in 1:20){
    #print(paste("now ", i, " is the optimal amino acid"))
    res[[i]] = eigOp(i,s,DisMat,MuMat,bf,C,Phi,q,Ne)
  }
  return(as.matrix(res))
}

ll3 <- function (dat1, tree, bf = c(0.25, 0.25, 0.25, 0.25), g = 1,
                 Q = c(1, 1, 1, 1, 1, 1), eig = NULL, assign.dat = FALSE,
                 ...)
{
  if (is.null(attr(tree, "order")) || attr(tree, "order") ==
    "cladewise")
    tree <- reorderPruning(tree)
  q = length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m = length(edge) + 1
  dat = vector(mode = "list", length = m)
  dat[1:q] = dat1[tree$tip.label]
  if (is.null(eig))
    eig = edQt(bf = bf, Q = Q)
  el <- tree$edge.length
  P <- getP(el, eig, g)
  nr <- as.integer(attr(dat1, "nr"))
  nc <- as.integer(attr(dat1, "nc"))
  node = as.integer(node - min(node))
  edge = as.integer(edge - 1)
  nTips = as.integer(length(tree$tip))
  mNodes = as.integer(max(node) + 1)
  contrast = attr(dat1, "contrast")
  nco = as.integer(dim(contrast)[1])
  res <- .Call("LogLik4", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
  result = res[[2]][[1]] + log(res[[1]][[1]] %*% bf)
  if (assign.dat) {
    dat[(q + 1):m] <- res
    attr(dat, "names") = c(tree$tip.label, as.character((q +
      1):m))
    assign("asdf", dat, envir = parent.frame(n = 1))
  }
  result
}
## For all the distinct sites (m), find loglikelihood, result is m * 20 matrix
## each column stores the loglikelihods when the corresponding amino acid is optimal
llaa <- function(tree,data,eigAll,bf=rep(1/20,20),C=2,Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  for(i in 1:20){ #when optimal aa is i
    llopi = ll3(data,tree,bf=bf,g=1,eig=eigAll[[i]])
    result = cbind(result,llopi)
  }
  opaa = apply(result,1,which.max)
  sitelik = apply(result,1,max)
  loglik = sum(weight * sitelik)
  return(list(loglik= loglik, llmat=result))
}

## find the loglikelihood given Q and other paramters, here Q is the lower triangular part of the 
## nucleotide transition rate matrix of length 6
mll <- function(data,tree,s,beta,gamma,Q=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  call <- match.call()
  if(class(tree)!="phylo") stop("tree must be of class phylo") 
  if (is.null(attr(tree, "order")) || attr(tree, "order") == 
    "cladewise") 
    tree <- reorderPruning(tree)
  if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
  if(is.null(tree$edge.length)) stop("tree must have edge weights") 
  if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")
  if(is.null(Q)) Q = rep(1,6)
  if(is.null(bfaa)) bfaa = findBf2(data)
  dismat = GM_cpv(GM_CPV,al,beta,gamma)
  mumat = aa_MuMat_form(Q)
  mumat = sym.to.Q(mumat,bfaa)
  eigall = eigAllaa(s,dismat,mumat,bf=bfaa,C,Phi,q,Ne)
  ll = llaa(tree,data,eigall,bfaa,C,Phi,q,Ne)
  result = list(ll=ll,data=data,tree=tree,s=s,GMweights=c(al,beta,gamma),Q=Q,bfaa=bfaa,call=call)
  class(result) = "mllm"
  return(result)
}