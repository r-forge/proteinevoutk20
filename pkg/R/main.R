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
NU_VEC <- c(1.96575,4.08655,1.39431,1.46172,6.36024,1.00000)
NU_VEC_JC <- rep(1.0,6)
##empirical base frequencies
BF <- c(0.31065,0.18365,0.20955,0.29615)
##Best tree with branch lengths output by PAUP
ROKAS_TREE <- read.nexus(paste(datadir,"GTR.tre",sep=""))
##Read the properties of amino acids, including c(composition),p(polarity) and 
##v(molecular volume), save the data.frame
GM_CPV <- read.csv(paste(datadir,"Grantham_cpv.csv",sep=""),header=TRUE, sep=",",row.names=1)
GM_CPV <- data.matrix(GM_CPV)

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
al <- avg_dis(GM_CPV[,1])
be <- avg_dis(GM_CPV[,2])
ga <- avg_dis(GM_CPV[,3])


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
GM <- GM_cpv(GM_CPV,al,be,ga)


##################################################################
##   construct matrices for mutation and substitution  ##
##################################################################
##Short names for all the amino acids, in alphabetical order
AA <- s2c("ARNDCQEGHILKMFPSTWYV") #amino aicds in alphabetical order
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
CDS <- nu_list[!nu_list %in% stop_cd] #61 non-stop codons
rm(nu_list)
rm(stop_cd)

##################################################################
##   Read gene data ##
##################################################################
## read in fasta file of nucleotide, convert it into amino acid data
## type = "num" -- integers ; "AA" -- amino acid names; "phyDat" --- phyDat data type used in phangorn
conv <- function(filename,type="num"){
  levels <- AA #amino acids in alphabeticla order (not the single letter names, the 3-letter names)
  data <- read.fasta(file=filename) #read fasta file including nucleotide data
  dna.to.aa <- function(x){ #given the order of the data, convert nucleotide data to amino acid data
    aas <- translate(seq=data[[x]])
    if(type=="num"){
      aas <- as.numeric(factor(aas,levels)) #convert from letters to numbers, as levels
      return(aas)
    }
    return(aas)
  }
  seqdata <- t(sapply(1:length(data),dna.to.aa,simplify="array")) #sapply and simplify to array
  dimnames(seqdata)[[1]]<- attr(data,"name") #Change the dimnames to the species names
  if(type=="phyDat") ## convert amino acid sequences to phyDat type
    seqdata = phyDat(seqdata,type="AA")
  return(seqdata)
}

##   Read in 106 gene data  from Rokas's ##
# l <- 106
# ROKAS_DATA <- vector("list",length=l)
# ##data stores the data of all genes in rokas's data
#   for(i in 1:l){
#     ROKAS_DATA[[i]] <- conv(paste(datadir,"gene",i,".fasta",sep=""))
#   }
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
## Find the most frequent element of a vector
## It works for numbers and characters
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#############################################################################
##given lower triangular part of R (i.e. Q = R %*% diag(bf)), and base frequencies, find the scaled Q
## Default: 4 by 4 matrix for Jukes-Cantor model
mat_form_lowtriQ <- function(Q=rep(1,6),bf=rep(0.25,4)){
  l=length(bf)
  res = matrix(0, l, l)
  res[lower.tri(res)] = Q
  res = res+t(res)
  res = res * rep(bf,each=l)
  diag(res) = -rowSums(res)
  res2 = res * rep(bf,l)
  diag(res2)=0 
  res = res/sum(res2)
  return(res)
}
## given Q and bf, find scaled Q
## Goal: sum(pi_i * Q_ii) = -1
scaleQ <- function(Q,bf){
  scalefactor = - sum(diag(Q)*bf) 
  return(Q/scalefactor)
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
codon_diff <- function(cd1,cd2){
  cd1 = s2c(cd1)
  cd2 = s2c(cd2)
  diff_ind = (cd1 != cd2)
  if(sum(diff_ind)==1)
    return(c(cd1[diff_ind],cd2[diff_ind]))
  else
    return(NULL)
}
#Find the mutation rate matrix (20 by 20) from the 6 rates in GTR mutation rate
#matrix for nucleotides
#the rate G<->C is normalized to 1, now the scale is to scale the rates to real values
aa_MuMat_form <- function(vec=rep(1,6),bf=rep(0.25,4)){
  Q <- mat_form_lowtriQ(vec,bf) #rate matrix for nucleotides    
  dimnames(Q) = list(Nu,Nu)
  ##mutation matrix for 61 codons
  codon_array <- array(0,dim=c(61,61),dimnames=list(CDS,CDS))
  for(m in 1:61){                       #loop through all 61 codons
    cd_str <- CDS[m] #mth codon in string format
    cd <- s2c(cd_str)   #one codon, convert from string to char vector
    for(i in 1:3){           #every nucleotide in the codon can mutate
      ngb <- cd #neighboring codon, set it equal to the starting codon for now
      for(j in 1:4){
        if(Nu[j]!=cd[i]){ #nucleotide has to change to a different one
          ngb[i] <- Nu[j] # change the ith position in codon to jth nucleotide
          a_new <- translate(ngb)       #new amino acid
          ngb_str <- c2s(ngb)
          if(a_new != "*")
            codon_array[cd_str,ngb_str] <- codon_array[cd_str,ngb_str] + Q[cd[i],Nu[j]]
        }
      }
    }
  }
  diag(codon_array) <- -rowSums(codon_array) #row sums equal to 0
  ##20*20 array: mutations rates between amino acids
  arr <- array(0,dim=c(20,20),dimnames=list(AA,AA)) #mutation rate matrix for amino acids
  for(i in 1:61){
    for(j in 1:61){
      from_codon <- CDS[i]
      to_codon <- CDS[j]
      from_aa <- translate(s2c(from_codon))
      to_aa <- translate(s2c(to_codon))
      arr[from_aa,to_aa] <- arr[from_aa,to_aa] + codon_array[from_codon,to_codon]
    }
  }
  #dimnames(arr) <- NULL
  return(arr)
}
## Amino acid mutation rate matrix for the GTR model of nucleotide
#MUMAT <- aa_MuMat_form(NU_VEC,BF)
#MUMAT_JC <- aa_MuMat_form()
#############################################################################
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
      if(MuMat[i,j]!=0){
        mat[i,j] <- 2*Ne*MuMat[i,j]*fix_protein(i,j,aa_op,s,DisMat,C,Phi,q,Ne) #fixation prob -> transition rate
        mat[j,i] <- 2*Ne*MuMat[j,i]*fix_protein(j,i,aa_op,s,DisMat,C,Phi,q,Ne) #symmetric entry
      }
    }#end for j
  }#end for i
  diag(mat) <- -rowSums(mat)
  return(mat)
}
######################################################
##likelihood function for 1 site
## Q: normalized (scaled) rate matrix
ll_site <- function(tree,data,Q,bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  ##If the given tree is not rooted and binary, then throw error and exit
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  m = 20
  ##if the base frequencies are not specified, do a uniform distribution
  if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
  GM1 = GM_cpv(GMcpv,alpha,beta,gamma)
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