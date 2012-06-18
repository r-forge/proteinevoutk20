##load the required libraries
library(seqinr)
library(phangorn)
library(optimx)
library(expm)
library(parallel)
library(minqa)
##all the names of the nucleotides and codons that'll be needed 

##Nucleotide mutation rates estimated from Rokas's data (using PAUP):
##Model GTR, with Gamma rate distribution and invariable sites and etc
##Nuvec <-  c(2.94194,8.23705,1.39133,2.33991,14.86435,1.00000)
##Model GTR
Nu_vec <- c(1.96575,4.08655,1.39431,1.46172,6.36024,1.00000)
##Nuvec <-  c(2.94194,8.23705,1.39133,2.33991,14.86435,1.00000)
freq <- c(0.31065,0.18365,0.20955,0.29615)

##Short names for all the amino acids, in alphabetical order
aa <- s2c("SRLPTAVGIFYCHQNKDEMW")
aa <- levels(factor(aa)) #all the amino acids
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

## read in fasta file of nucleotide, convert it into amino acid data
conv <- function(filename){
  levels <- levels(factor(s2c("SRLPTAVGIFYCHQNKDEMW"))) #amino acids in alphabeticla order
  data <- read.fasta(file=filename) #read fasta file including nucleotide data
  seqdata <- vector()
  for(i in 1:length(data)){
    aas <- translate(seq=data[[i]]) #translate to amino acids
    aan <- as.numeric(factor(aas,levels)) #convert them to numbers
    seqdata <- rbind(seqdata,aan) #add the amino acids for this species to the result vector
  }
  dimnames(seqdata)[[1]]<- attr(data,"name") #Change the dimnames to the species names
  return(seqdata)
}

##rokas.gene1 <- conv("RokasGene1.fasta") #Read in the data of 1st gene in Rokas's data
##rokas.gene1 <- conv("gene1.fasta") #Read in the data of 1st gene in Rokas's data

##Given sequence data on a gene for different species, resample data of the 
##same length for bootstrapping 
##PAY SPECIAL ATTENTION TO THE SEED OF RANDOM GENERATOR FUNCTION!!!
##LOOK AT THE BOOT PACKAGE AND SEE IF THERE IS SOME ALREADY EXISTED FUNCTION
resample <- function(seqdata){
  l <- dim(seqdata)[2] #number of (amino acid) sites in the gene
  sample_index <- sample(l,l,replace=T) #sample the indices with replacements
  SeqSample <- seqdata[,sample_index] #return the new sample
  dimnames(SeqSample) = dimnames(seqdata)
  return(list(SeqSample=SeqSample, index=sample_index))
}

## Find the most frequent element of a vector
## It works for numbers and characters
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
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
  rate_mat[lower.tri(rate_mat)] <- rate_mat[upper.tri(rate_mat)]
  for(i in 1:4){ #row sums equal to 0
    rate_mat[i,i] <- -sum(rate_mat[i,-i])
  }
  dimnames(rate_mat)[[1]] <- Nu
  dimnames(rate_mat)[[2]] <- Nu
  return(rate_mat)
}

##Find the scale for the subsitution rate matrix Q
scale.rate <- function(Q){
  Pi <- as.vector(expm(Q*10000)[1,]) #equilibrium frequency vector
  return(-sum(Pi*diag(Q)))
}
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
  for(i in 1:61){                       #row sums equal to 0
    codon_array[i,i] <- -sum(codon_array[i,-i])
  }
  #browser()
  ##equilibrium frequenciew for 61 codons
  cd_eqlm <- expm(codon_array*10000)[1,]
	
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
  return(arr)
}

##vec <- runif(6)
##aa_MuMat_form(vec)


##setwd("~/Documents/ProteinEvo") #set working directory
##Read the properties of amino acids, including c(composition),p(polarity) and 
##v(molecular volume), save the data.frame
GMcpv <- read.csv("Grantham_cpv.csv",header=TRUE, sep=",",row.names=1)
GMcpv <- data.matrix(GMcpv)
dimnames(GMcpv) <- NULL
##Build the Gramtham matrix with coefficients for all components as parameters
##alpha for composition,beta for polarity,gamma for molecular volume
##The mean of the distance is normalized to 1,(mean of only the upper triangle elements of the matrix)
##The parameter values in Grantham matrix are:
##alpha=1.833, beta=0.1018,gamma=0.000399

##this function computes the weights that are used by Grantham
## (dispite the name, they are not the average distances though, instead they are functions of average distances)
avg_dis <- function(vector){
  l <- length(vector)
  Sum <- 0
  for(i in 1:(l-1)){
    Sum <- Sum + sum(abs(vector[i]-vector[(i+1):l]))
  }
  ct <- l*(l-1)/2
  return((ct/Sum)^2)
}
##Grantham weights
al <- avg_dis(GMcpv[,1])
be <- avg_dis(GMcpv[,2])
ga <- avg_dis(GMcpv[,3])

#Given the data from Grantham on all three factors, and the weights,
#Find the distance matrix for all the amino acids (20 by 20) with mean 1
#Grantham matrix has mean 100 instead of 1
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
##GM <- GM_cpv(GMcpv)


##################################
##physiochemical distance matrix#
##################################
##parameter values, now passed into the functions instead of being global
## a1 <- 2
## a2 <- 1
## Phi <- 0.5
## q <- 4*10^(-3)
## Ne <- 1.37*10^3


#############################################################################
##only the first m states available for each site, instead of all 20###
## now m is passed into functions as parameters instead of global variables
#############################################################################
##find the physiochemical distance vector between two proteins, given the distance matrix
pchem_d <- function(protein1, protein2,DisMat){
  l <- length(protein1)
  d <- vector(length=l)
  site_d <- function(k){
    i <- protein1[k]
    j <- protein2[k]
    return(DisMat[i,j])
  }
  d <- sapply(c(1:l),site_d,simplify=TRUE) 
  #with sapply, return a vector or array by default
  #with apply, return a list by default
  return(d)
}

#################################################################
#Functionality function, arguments: 
#d: vector of distances between the optimal and observed amino acids
#s: seleciton coeffecient, can be passed as a vector or a number if it's the same for all sites
Ftny <- function(d, s, model=2){  
  if((length(s)==1)&&(length(d)!=1)){ #if s is given as a scalar, then treat it to be the same across all sites
    s <- rep(s,length(d))
  }
  else{
  }
  if(model==1){ #model 1
    result <- prod(exp(-d*s)) #product of ftny at sites (exponential functions)
  }
  else{ #model 2
    result <- length(d)/(sum(1+d*s)) #harmonic mean of ftny at all sites (F_i = 1 + d_i * s_i)
  } 
  return(result)
} #d and s are given

#Functionality of a protein given the optimal protein (aa sequence) and parameters
Ftny_protein <- function(protein,protein_op,s,DisMat,model=2){
  d <- pchem_d(protein,protein_op,DisMat)
  return(Ftny(d,s,model))
}
#fitness function
#d1, d2: distance vectors
#s: selection (vector or scalar)
#Ne: effective population size
#q, Phi: constants
#C: cost function, linearly function of the length of the protein
#C_n = a1 + a2*n, cost --  linear function of the length
#formula number: 1: Sela & Hirsh's formula (default); 2: canonical formula
fix <- function(d1, d2, s, indep=TRUE,formula=1, model=2,a1=2, a2=1, Phi=0.5,q=4e-7,Ne=1.36e7){
  if(length(d1)!=length(d2))
    return("error: two proteins have different lengths")
  n <- length(d1)
  C_n <- a1 + a2*n
  if(n==1){ #if there is only 1 amino acid in the sequence
    if((d1==d2)||(s==0)){ #When the fitnesses are the same, neutral case, pure drift
      return(1/(2*Ne))
    }
    else{
      ftny1 <- Ftny(d1,s,model)
      ftny2 <- Ftny(d2,s,model)
      fit_ratio <- exp(-q*Phi*C_n*((1/ftny1)-(1/ftny2))) #w1/w2
    }
  }
  else{ # there are more than 1 amino acids in the protein
    if((length(d1)!=1)&&(length(s)==1)){
      s <- rep(s,n)
    }
    index <- c(1:n)
    diff <- d1 - d2
    S <- index[diff==0]
    D <- index[diff!=0]
    if(length(D) > 1){ #if there are more than 1 site with different amino acids
      return(0)
    }
    else if(length(D)==0){ #two protein sequences are the same
      return(1/(2*Ne)) 
    }
    else{ #if there is only one site that is different
      if(model==1){
        ftny_D1 <- Ftny(d1[D],s[D],model) #different terms in the product
        ftny_D2 <- Ftny(d2[D],s[D],model)
        if(ftny_D1==ftny_D2){
          return(1/(2*Ne))
        }
        else{
		if(indep==TRUE){ # in site-indepdent case, F_S is equivalent to 0
          		ftny_S <- 1
        	}
        	else{
          		ftny_S <- Ftny(d1[S],s[S],model) #same terms in the product for both proteins
        	}
          fit_ratio <- exp(-q*Phi*C_n*(1/ftny_S)*((1/ftny_D1)-(1/ftny_D2)))
        }
      } #end if model==1
      else{ #model==2, results are the same for dependent and independent cases
		#now the length of D should be equal to 1
        if(s[D]==0)
          return(1/(2*Ne))
        else
          fit_ratio <- exp(-q*Phi*C_n*(d1[D]-d2[D])*s[D]/n)
      } #end if model==2
    }
  }
  #different fixation prob formula
  if(fit_ratio==1){
    return(1/(2*Ne))
  }
  else if(fit_ratio==Inf){ #Is this true?
    return(0)
  }
  else{
      if(formula==1){
        result <- (1-fit_ratio)/(1-fit_ratio^(2*Ne)) #diploid in W-F process
      }
      else if(formula==2){
        se <- 1/fit_ratio - 1 # s=(w2-w1)/w1, selection advantage of 2 comparing to 1
        #result <- (1-exp(-2*se))/(1-exp(-4*Ne*se))
        result <- (1-exp(-se))/(1-exp(-2*Ne*se))
      }
  }
  return(result)
}

#fixation probability of a mutant with initial freq 1/(2Ne) in a diploid population
#with proteins as parameteres instead of distance vectors for the proteins
fix_protein <- function(protein1, protein2, protein_op, s, DisMat,
                        indep=FALSE, formula=1, model=2,a1=2, a2=1, Phi=0.5,q=4e-7,Ne=1.36e7){
  d1 <- pchem_d(protein1,protein_op,DisMat)
  d2 <- pchem_d(protein2,protein_op,DisMat)
  return(fix(d1,d2,s,indep,formula,model,a1,a2,Phi,q,Ne))
}

######################################################

#Function to generate the transition matrix (m*m for now) 
#given the optimal amino acid. This is only for one site (amino acid)
#parameter: aa_op: optimal amino acid at this site
#s: selection coefficient(SCALAR, since it's only one site)
#mu: mutation rate from one aa to another, they are assumed to be the same across all amino acids.
#m: number of states (amino acids) considered
#DisMat: distance matrix between amino acids
#MuMat: mutation rate matrix between amino acids
mat_gen_indep <- function(aa_op,m,s,DisMat,MuMat,formula=1,model=2,a1=2, a2=1, Phi=0.5,q=4e-7,Ne=1.36e7){
  mat <- matrix(0,nrow=m,ncol=m) #set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      mat[i,j] <- 2*Ne*MuMat[i,j]*fix_protein(i,j,aa_op,s,DisMat,indep=T,formula,model,a1,a2,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- 2*Ne*MuMat[i,j]*fix_protein(j,i,aa_op,s,DisMat,indep=T,formula,model,a1,a2,Phi,q,Ne) #symmetric entry
    }#end for j
    mat[i,i] <- -(sum(mat[i,])) #diagonal entry; sum of each row is 0
  }#end for i
  mat[m,m] <- -sum(mat[m,])
  return(mat)
}



##Given a protein, the optimal protein, and the selection vector (or constant selection)
##find the vector of rates of moving from the protein to all its neighbors
##If the length of the protein is l, then it has (m-1)*l neighbors
rate_move <- function(protein, protein_op, m, s,DisMat,MuMat,indep=TRUE,formula=1,model=2,a1=2, a2=1, Phi=0.5,q=4e-7, Ne=1.36e7){
  l <- length(protein) #number of sites
  rates <- vector(length=l*(m-1)) # rates of moving to neighboring proteins, vector to return
  for(i in 1:l){
	aa <- c(1:m)[-protein[i]]
	mu_vec <- MuMat[protein[i],aa]
	rates[((i-1)*(m-1)+1):(i*(m-1))] <- 2*Ne*mu_vec*mapply(fix_protein,rep(protein[i],(m-1)),aa,rep(protein_op[i],m-1),MoreArgs=list(s=s,DisMat=DisMat,indep=indep,
				formula=formula,model=model,a1=a1,a2=a2,Phi=Phi,q=q,Ne=Ne),SIMPLIFY=TRUE)
  }
  rates
}



#Simulation. Given the starting protein and the optimal protein
#t: running time of the chain
#MuMat: matrix of mutation rates between amino acids
simulation <- function(protein,protein_op,t,m,s,DisMat,MuMat,indep=TRUE,formula=1,model=2,a1=2, a2=1, Phi=0.5,q=4e-7, Ne=1.36e7,record.fitness=FALSE){
  #browser()
  l <- length(protein) #number of sites
  C <- a1 + a2*l #protein production cost
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #colomn names
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  while(t_now < t){ #when current time is less than t
    rates <- rate_move(protein,protein_op,m,s,DisMat,MuMat,indep,formula,model,a1,a2,Phi,q,Ne) #moving rates of a protein to its neighbors
    nonzero_index <- which(rates!=0) #indices where the rates are not equal to 0
    rates <- rates[nonzero_index] #only keep the nonzero entries
    lambda <- sum(rates) # rate of staying at the same state
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait
    if(t_now < t){
      #rates <- rates/lambda #probability vector of moving
      r <- runif(1)
      index <- sample(nonzero_index,1,replace=T,rates) #index of the protein it moves to
      pos <- (index-1) %/% (m-1) + 1 #index of the site that changes
      rmd <- (index-1) %% (m-1) + 1 # the amino acid that site changes to
      protein[pos] <- seq(1:m)[-protein[pos]][rmd] #new protein
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  if(record.fitness){
    ftny_vec <- NULL
    fitness_vec <- NULL
    for(i in 1:dim(path)[1]){
      ftny_vec <- c(ftny_vec,Ftny_protein(path[i,seq(1:l)],protein_op,s,DisMat,model)) #functionality
    }
    fitness_vec <- exp(-q*Phi*C/ftny_vec) #fitness
    path <- cbind(path,ftny_vec)
    path <- cbind(path,fitness_vec)
    colnames(path)[l+3] <- "Functionality"
    colnames(path)[l+4] <- "Fitness"
  }
  ##shift the third column up one step so that the waiting time is the time
  ## spent in the state in the first column
  path[,3] <- c(path[-1,3],NA)  
  path
}
#simulation(c(1,2),c(3,4),1000,10,0.1,GM,mumat,indep=T)


#Simulation of protein sequences of length "l" on a phylogeny "tree", 
#given the ancestral sequence "rootseq", optimal amino acid sequence "protein_op",
#selection coefficient "s", 
simTree <- function(tree,l,protein_op=rep(1,l),m,s,GTRvec,alpha=al,
                    beta=be, gamma=ga,scale=T,indep=TRUE,formula=1,model=2,a1=2, a2=1, 
                    Phi=0.5,q=4e-7, Ne=1.36e7,
                    bf=NULL,rootseq=NULL,ancestral=FALSE){
    # check if the input tree is rooted binary tree. if not, throw an error.
    if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
    if(is.null(bf)) bf = rep(1/m,m) #base frequency, randomly chosen from all states
    if(is.null(rootseq))rootseq = sample(c(1:m), l, replace=TRUE, prob=bf) #sequence at the root
    GM1=GM_cpv(GMcpv,alpha,beta,gamma) #Distance matrix from new weights
    mumat = aa_MuMat_form(GTRvec,scale)
    tree = ape:::reorder.phylo(tree) #oder the tree, cladwise
    edge = tree$edge #edges
    nNodes = max(edge) #number of nodes in the tree
    res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes
    parent <- as.integer(edge[, 1]) #parents of the edges
    child <- as.integer(edge[, 2]) #children of the edges
    root <- as.integer(parent[!match(parent, child, 0)][1]) #root node
    res[root,] <- rootseq
    tl = tree$edge.length #lengths of the edges
    #browser()
    for(i in 1:length(tl)){
        from = parent[i] 
        to = child[i]
    ##simulation on this one branch
        for(j in 1:l){
        res[to,j] = as.numeric(tail(simulation(res[from,j],protein_op[j],tl[i],m,s,GM1,mumat,indep,
                              formula,model,a1,a2,Phi,q,Ne,record.fitness=FALSE),1)[1])
        
        }
    }
    k = length(tree$tip)
    label = c(tree$tip, as.character((k+1):nNodes))
    rownames(res)=label 
    if(!ancestral)res = res[ tree$tip, , drop=FALSE]
    return(as.data.frame(res))
    #if(pt=="AA") return(phyDat.AA(as.data.frame(res), return.index=TRUE))
}
#simTree(tree,5,protein_op=Protein_op[1:5],m=20,s=0.3,Nu_vec,al,be,ga,q=4e-7,Ne=1.36e7,Root[1:5])
########################################################
#Find the likelihood of a tree given data
########################################################
#data with only 1 site

##A function I need to do the matrix exponential correctly
is.nan.inf <- function(x){
  return((any(x==Inf))||(any(is.nan(x))))
}
## A modified version of expm that make sure there is no Inf and NaN in the result of matrix exponential
expm.m <- function(x){
  edefault <- expm(x) #expm from expm package with default method "higham08.b"
  if(is.nan.inf(edefault)==FALSE) #if there is no Inf or NaN values in the result, return the result
    return(edefault)
  else{
    eward77 <- expm(x,method="Ward77") #method Ward77
    if(is.nan.inf(eward77)==FALSE)
      return(eward77)
    else{
      epade <- expm(x,method="Pade") #method Pade
      if(is.nan.inf(epade)==FALSE)
        return(epade)
      else{
        etaylor <- expm(x,method="taylor") #method Taylor
        if(is.nan.inf(etaylor)==FALSE)
          return(etaylor)
        else{
          return(edefault)
          warning("all methods give Inf or NaN")
        }
      }
    }
  }
}

ll_site <- function(tree,data,optimal,m=20,s=1,MuMat,alpha=al, beta=be, gamma=ga,
                    root=NULL,bf=NULL,formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.36e7){
    ##If the given tree is not rooted and binary, then throw error and exit
    if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
    ##If the amino acid at the root is specified, set the base frequencies
    if(!is.null(root)){ 
      bf = rep(0,m)
      bf[root]=1 #set the frequency of the root aa to be 1, others to be 0
    }
    ##if the base frequencies are not specified, do a uniform distribution
    if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
    GM1 = GM_cpv(GMcpv,alpha,beta,gamma)
    ##MuMat = aa_MuMat_form(GTRvec)
    Q = mat_gen_indep(optimal,m,s,GM1,MuMat,formula,model,a1,a2,Phi,q,Ne) #transition rate matrix for the site, given the optimal aa
    
    tree <- ape:::reorder.phylo(tree,"p") #reorder the tree in pruningwise order
    edge = tree$edge #edges
    nNodes = max(edge) #number of nodes in the tree
    probvec = matrix(NA,nNodes,m) #probability of gettting different states at nodes that evolve to the current sequences
    parent <- as.integer(edge[, 1]) #parents of the edges
    child <- as.integer(edge[, 2]) #children of the edges
    root <- as.integer(parent[!match(parent, child, 0)][1])  
    tip <- as.integer(child[!match(child,parent,0)])
    for(i in 1:length(tip)){
      probvec[i,] = 0
      probvec[i,data[i]] = 1
    }
    #or use this to find tips:
    #tip <- c(1:length(tree$tip.label))
    tl = tree$edge.length #lengths of the edges
    for(i in 1:tree$Nnode){
        from = parent[2*i] 
        to = child[(2*i-1):(2*i)]
        t_left = tl[2*i-1]
        t_right = tl[2*i]
        ##pick the right method to calculate the matrix exponential
        v.left <- expm.m(Q*t_left)
        v.right <- expm.m(Q*t_right)
        for(j in 1:m){
          #find the likelihood of a node when likelihood of both its descendents is found --- pruning
          probvec[from,j] = (v.left[j,] %*% probvec[to[1],])*(v.right[j,] %*% probvec[to[2],]) 
        }
        #browser()
        if(sum(probvec[from,])==0) #probability is very very low
          warning("numerical overflow")
    }
    return(as.numeric(probvec[root,] %*% bf))
}
#ll_site(t,data[,1],optimal[1],20,0.1,aamat)
#ll_site(tree,data[[12]][,4],Mode(data[[12]][,4]),m=20,s=70,mumat,al,be,ga)

##Likelihood of data on a tree, given the selection coefficient "s" and the optimal amino acid sequence "protein_op"
##If protein_op (optimal protein) is not given, get it from the most frequent amino acids appeared in the data
ll_indep <- function(s,alpha,beta,gamma,MuMat,tree,data,m=20,protein_op=NULL,root=NULL,bf=NULL,
                             formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.36e7){
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  l <- dim(data)[2] #number of sites in the sequence
  pr <- 0
  ##The following only works with trees of 8 tips, revise accordingly, or change to a better way!!!!
  data <- data[,order(data[1,],data[2,],data[3,],data[4,],data[5,],data[6,],data[7,],data[8,])] #order the data so that same data will be next to each other
  if(is.null(protein_op)) #Use the most frequent amino acid, what if there are more than 1 with highest frequency??
    protein_op <- apply(data,2,Mode)
  llh = log(ll_site(tree,data[,1],protein_op[1],m,s,MuMat,alpha,beta,gamma,root,bf,formula,model,a1,a2,Phi,q,Ne))
  #print(llh)
  for(i in 2:l){
    if(length(which(data[,i]!=data[,i-1]))==0){ #if data at the two sites are the same, don't need to calculate again
      pr = pr - llh
      #print(llh)
    }
    else{ #if this pattern hasn't appeared before, calculate the loglikelihood
    llh = log(ll_site(tree,data[,i],protein_op[i],m,s,MuMat,alpha,beta,gamma,root,bf,formula,model,a1,a2,Phi,q,Ne))
    pr = pr - llh
    #print(llh)
    }
  }
  pr
}

##This function finds the MLE estimator for "s" only, given that all other parameters are known,
##including weights(alpha, beta, gamma), mutation rates, C, q, Phi, and Ne
MLE_GTR <- function(start_pt,lowerb, upperb,tree,data,m=20,alpha,beta,gamma,MuMat,optim.m=1,
                    protein_op=NULL,root=NULL,bf=NULL,formula=1,model=2,
                       a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.37e7){
  negloglike <- function(s){ #The function to minimize
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,formula,model,a1,
                    a2,Phi,q,Ne))
  }
  #ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,all.methods=T,save.failures=T,
                #hessian=T,control=list(trace=1))
  if(optim.m==1)
    ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,method="nlminb",
                hessian=F,control=list(trace=1))
  else #Powell method bobyqa
    #ans <- bobyqa(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(iprint=3))
    ans <- bobyqa(start_pt,negloglike,lower=lowerb,upper=upperb)#does not print out too much information
  #ans <- nlminb(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(trace=1))
  return(ans)
}
##Use the GS (golden section) method for optimization since there is only one parameter to optimize
MLE_GTR_GS <- function(start_pt,lowerb, upperb,tree,data,m=20,alpha,beta,gamma,MuMat,
                    protein_op=NULL,root=NULL,bf=NULL,formula=1,model=2,
                       a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.37e7){
  negloglike <- function(s){ #The function to minimize
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,formula,model,a1,
                    a2,Phi,q,Ne))
  }
  ans <- GS(negloglike,lowerb,start_pt,upperb)
  return(ans)
}
##MLE_GTR(tree,data,rep(2,10),10,0.1,al,be,ga)

##For given beta and gamma, (suppose that alpha=al or 1) find MLE for s in each gene
## and the maximum likelihood values under these s values
ll_sw <- function(beta, gamma, tree,data,alpha,MuMat,protein_op=NULL,root=NULL,bf=NULL,
                  formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.37e7,m=20){
  ##Do this for every gene in the data
  MLE_GTR(0.3,-10,1000,tree,data,m,alpha,beta, gamma,MuMat,protein_op,root,bf,formula,
          model,a1,a2,Phi,q,Ne)
}
##This function finds the MLE estimators for s, and the cpv weights. Here alpha is given, beta and gamma are being
##estimated, since only 2 of them are free parameters.
##Brian suggested that the restriction should be that alpha + beta + gamma = 1, DO THIS LATER!
MLE_GTR_sw<- function(start_pt,lowerb,upperb,tree,data,m=20,alpha,MuMat,
                      protein_op=NULL,root=NULL,bf=NULL,formula=1,model=2,
                      a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.37e07){
  negloglike <- function(para){
    s = para[1]
    beta = para[2]
    gamma = para[3]
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,root,bf,formula,model,a1,
                    a2,Phi,q,Ne))
  }
  #ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,method="nlminb",hessian=T,control=list(trace=1))
  ans <- bobyqa(start_pt,negloglike,lower=lowerb,upper=upperb,control=list(iprint=3))
  return(ans)
}


##This function finds the MLE's for s, two of the cpv weights, and the effective population size Ne
##Since the scale of Ne is big, we try to find the MLE for log_10(Ne)
MLE_GTR_swNe<- function(start_pt,lowerb,upperb,tree,data,m=20,alpha,MuMat,protein_op=NULL,
                        root=NULL,bf=NULL,formula=1,model=2,
                        a1=2,a2=1,Phi=0.5,q=4e-3){
  
  negloglike <- function(para){
    s = para[1]
    beta = para[2]
    gamma = para[3]
    Ne = 10^(para[4])
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,protein_op,
                    root,bf,formula,model,a1,a2,Phi,q,Ne))
  }
  browser()
  ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,method="nlminb",hessian=T,
                control=list(trace=1))
  return(ans)
}

##MLE's of two cpv weights -- beta and gamma
MLEweights <- function(start_pt,lowerb,upperb,tree,data,m=20,s,alpha,protein_op=NULL,root=NULL,bf=NULL,
                       formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-3,Ne=1.37e03){
  negloglike <- function(cpv){
    beta <- cpv[1]
    gamma <- cpv[2]
    return(ll_indep(s,alpha,beta,gamma,MuMat,tree,data,m,
                    protein_op,root,bf,formula,model,
                    a1,a2,Phi,q,Ne))
  }
  ans <- optimx(start_pt,negloglike,lower=lowerb,upper=upperb,method="nlminb",hessian=T,control=list(trace=1))
  return(ans)
}

##in this function, the root has equilibrium frequencies, which are calculated from
##Sella-Hirsh's formula, hence depend on s.
##In reality, the root doesn't affect the result a lot, given that the branches are long enough
##for all the sites to reach equilibrium
ll_root_eqlm <- function(s,alpha,beta,gamma,tree,data,protein_op,m=20,
                             formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-3,Ne=1.37e03,mu=1/(2*Ne)){
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  eqlm_freq_mat <- matrix(NA,m,m)
  C_n <- a1+a2
  for(k1 in 1:m){
    ftny_vec <- vector(length=m)
    for(k2 in 1:m){
      ftny_vec[k2] <- Ftny_protein(k2,k1,s,alpha,beta,gamma,model)
    }
    fit_vec <- exp(-C_n*Phi*q/ftny_vec)
    freq <- (fit_vec/max(fit_vec))^(2*Ne-1)
    eqlm_freq_mat[k1,] <- freq/sum(freq)
  }
  l <- length(protein_op)
  pr <- 0
  for(i in 1:l){
    llh <- ll_site(tree,data[,i],protein_op[i],m,s,alpha,beta,gamma,eqlm_freq_mat[protein_op[i],],formula,model,a1,a2,Phi,q,Ne,mu)
    pr <- pr - log(llh)
  }
  return(pr)
}

#find the MLE for s when the equilibrium frequencies are not given
#instead they come from Sella-Hirsh's formula and depend on s
MLEeqlm <- function(tree,data,protein_op,m=20,formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-7,Ne=1.37e07){
    negloglike <- function(s){
    return(ll_root_eqlm(s,tree,data,protein_op,m,formula,model,a1,a2,Phi,q,Ne))
  }
  best <- GS(negloglike,0,2)
  return(best)
}

GS <- function (func, a, b, tol=1e-5) 
{
  r <- (sqrt(5)-1)/2
  x <- b - r* (b-a)
  y <- a + r* (b-a)
  fx <- func(x)
  fy <- func(y)
  steps <- ceiling (log(2*tol/(b-a))/log(r))
  for (i in 1:steps) {
    if ((fy > fx)||(fy==Inf)) {
      b <- y
      y <- x
      fy <- fx
      x <- b - r* (b-a)
      print(paste("x=",x))
      fx <- func(x)
      print(paste(fx,", at",x))
    }
    else {
      a <- x
      x <- y
      fx <- fy
      y <- a + r* (b-a)
      print(paste("y=",y))
      fy <- func(y) 
      print(paste(fy,", at",y))
    }
  }
  
  if (fx < fy) {
    position <- x
    value <- fx
  }
  else {
    position <- y
    value <- fy
  }
  if ((abs((position-a))<tol)&&(func(position)>func(a))){
    position <- a
    value <- func(a)
  }
  else if((abs((position-b))<tol)&&(func(position)>func(b))){
    position <- b
    value <- func(b)
  }  
  
  response <- list(position=position, value=value) 
  return (response)
}

##Another version of golden section 
GS <- function(f, x1, x2, x3, tol=1e-7) {
  r <- 2 - (1+sqrt(5))/2
	x4 <- x2 + (x3-x2)*r
	if ( abs(x3-x1) < tol ){
		return(list(pos=x2,value=f(x2)))
	}
	if (f(x4) < f(x2)) {
		gSection(f, x2, x4, x3, tol)
	} else {
		gSection(f, x4, x2, x1, tol)
	}
}