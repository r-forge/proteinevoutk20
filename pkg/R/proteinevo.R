setwd("~/Documents/ProteinEvo") #set working directory
library("Matrix") #load the library that includes the function "expm" to calculate exponential of matrix
library("phangorn")


##################################
#physiochemical distance matrix#
##################################
#parameter values, now passed into the functions instead of being global
# a1 <- 2
# a2 <- 1
# Phi <- 0.5
# q <- 4*10^(-3)
# Ne <- 1.37*10^3
#############################################################################
#Import Grantham distance matrix from .csv file, and make the mean distance 
# equal to 1 instead of original 100
GranthamMatrix <- read.csv("GranthamMatrix.csv",
                           header=TRUE, sep=",",row.names=1)
#View(GranthamMatrix) #check the matrix
GM <- GranthamMatrix #store it in a different matrix and manipulate it
#complete the upper triangle matrix to whole, with diagonals equal to 0
GM[1,1] <- 0
for(i in 2:20){
  GM[i,i] <- 0
  for(j in 1:(i-1)){
    GM[i,j] <- GM[j,i] #symmetric matrix
  }
}
GM <- GM/100 #make the mean 1 instead of 100####

#############################################################################
#m <- 3 #only the first m states available for each site, instead of all 20###
# now m is passed into functions as parameters instead of global variables
#############################################################################
#find the physiochemical distance vector between two proteins
pchem_d <- function(protein1, protein2){
  l <- length(protein1)
  d <- vector(length=l)
  for(i in 1:l){
    d[i] <- GM[protein1[i],protein2[i]]
  }
  return(d)
}

#################################################################
#Functionality function, arguments: 
#d: vector of distances between the optimal and observed amino acids
#s: seleciton coeffecient, can be passed as a vector or a number if it's the same for all sites
Ftny <- function(d, s, model=1){  
  #browser()
  if(length(s)==1){ #if s is given as a scalar, then treat it to be the same across all sites
    s <- rep(s,length(d))
  }
  else{
  }
  if(model==1){ #model 1
    result <- prod(exp(-d*s)) #product of ftny at sites (exponential functions)
  }
  else if(model==2){ #model 2
    result <- length(d)/(sum(1+d*s)) #harmonic mean of ftny at all sites (F_i = 1 + d_i * s_i)
  } 
  #Other formulae?
  return(result)
} #d and s are given

Ftny_protein <- function(protein,protein_op,s,model=1){
  d <- pchem_d(protein,protein_op)
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
fix <- function(d1, d2, s, indep=FALSE,formula=1, model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3){
  if(length(d1)!=length(d2))
    return("error: two proteins have different lengths")
  n <- length(d1)
  C_n <- a1 + a2*n
  if(n==1){ #if there is only 1 amino acid in the sequence
    if((d1==d2)|(s==0)){ #When the fitnesses are the same, neutral case, pure drift
      return(1/(2*Ne))
    }
    else{
      ftny1 <- Ftny(d1,s,model)
      ftny2 <- Ftny(d2,s,model)
      fit_ratio <- exp(-q*Phi*C_n*((1/ftny1)-(1/ftny2))) #w1/w2
    }
  }
  else{ # there are more than 1 amino acids in the protein
    if(length(s)==1){
      s <- rep(s,n)
    }
    index <- seq(1,n)
    diff <- d1 - d2
    S <- index[diff==0]
    D <- index[diff!=0]
    if(length(D) > 1){ #if there are more than 1 site with different amino acids
      return(0)
    }
    else if(length(D)==0){ #two protein sequences are the same
      return(1/(2*Ne)) 
    }
    else if(length(D)==1){ #if there is only one site that is different
      if(model==1){
        if(indep==TRUE){ # in site-indepdent case, F_S is equivalent to 0
          ftny_S <- 1
        }
        else{
          ftny_S <- Ftny(d1[S],s[S],model) #same terms in the product for both proteins
        }
        ftny_D1 <- Ftny(d1[D],s[D],model) #different terms in the product
        ftny_D2 <- Ftny(d2[D],s[D],model)
        if(ftny_D1==ftny_D2){
          return(1/(2*Ne))
        }
        else{
          fit_ratio <- exp(-q*Phi*C_n*(1/ftny_S)*((1/ftny_D1)-(1/ftny_D2)))
        }
      } #end if model==1
      else{ #model==2, results are the same for dependent and independent cases
        ds_diff <- (d1[D]-d2[D])*s[D]
        if(ds_diff==0)
          return(1/(2*Ne))
        else
          fit_ratio <- exp(-q*Phi*C_n*ds_diff/n)
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
fix_protein <- function(protein1, protein2, protein_op, s, indep=FALSE, formula=1, model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3){
  d1 <- pchem_d(protein1,protein_op)
  d2 <- pchem_d(protein2,protein_op)
  return(fix(d1,d2,s,indep,formula,model,a1,a2,Phi,q,Ne))
}

######################################################

#Function to generate the transition matrix (m*m for now) 
#given the optimal amino acid. This is only for one site (amino acid)
#parameter: aa_op: optimal amino acid at this site
#s: selection coefficient(SCALAR, since it's only one site)
#mu: mutation rate from one aa to another, they are assumed to be the same across all amino acids.
#m: number of states (amino acids) considered
mat_gen_indep <- function(aa_op,m,s,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne)){
  mat <- matrix(0,nrow=m,ncol=m) #set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      mat[i,j] <- 2*Ne*mu*fix_protein(i,j,aa_op,s,formula,model,a1,a2,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- 2*Ne*mu*fix_protein(j,i,aa_op,s,formula,model,a1,a2,Phi,q,Ne) #symmetric entry
    }#end for j
    mat[i,i] <- -(sum(mat[i,])) #diagonal entry; sum of each row is 0
  }#end for i
  mat[m,m] <- -sum(mat[m,])
  return(mat)
}

#Given an array whose rows are proteins, find the substitution rate matrix between them
#arr: the protein sequences; each row represents a protein 
#protein_op: optimal protein
#s: selection coefficients 
#########################################################################3
mat_gen_arr<- function(arr,protein_op,s,indep=FALSE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3, mu=1/(2*Ne)){
  n1 <- dim(arr)[1] #number of proteins to consider
  #n <- dim(arr)[2] #number of sites in each protein
  mat <- matrix(0,nrow=n1,ncol=n1) #matrix to return, every entry is initiated to be 0
  for(i in 1:(n1-1)){
    for(j in (i+1):n1){
      mat[i,j] <- 2*Ne*mu*fix_protein(arr[i,],arr[j,],protein_op,s,indep,formula,model,a1,a2,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- 2*Ne*mu*fix_protein(arr[j,],arr[i,],protein_op,s,indep,formula,model,a1,a2,Phi,q,Ne) #symmetric entry
    }
    mat[i,i] <- -sum(mat[i,])
  }
  mat[n1,n1] <- -sum(mat[n1,])
  return(mat)
}
#######################################################################################################


#######################################################################################################
#Given 2 amino acid sequences, find the transition rate matrix using parsimony path
#1.Find all the sequences on the path
#2.Calculate the entries between 2 sequences that differ at only one site
#3.If more than 1 site differ, the rate is 0
#4.What about the optimal sequence? (Right now the optimal amino acid sequence is given
#If the number of sites that are differ is m, then first step, there are m choices to make, 
# second step, there are m-1 choices, then m-2 choices, hence the total number of pasimony paths
# is m!. How many sequences are possible at each step? 1st: m; 2nd: C(m,2); 3rd: C(m,3)..
# Hence in total the number of sequences need to consider is 2^m. Use binary numbers to do the
# count. From 0000..0 to 111..1
#Convert a decimal integer to a binary number
binary <- function(x,k){
  #ndigits <- (floor(logb(x, base=2))+1)
  Base <- rep(0,k)
  for(i in 1:k){#i <- 1
    Base[k-i+1] <- (x %% 2)
		x <- (x %/% 2)
	}
	return(Base)
}

integer.base.b <-function(x, b=2){
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
  	print(list(ERROR="x not integer", x=x))
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=b))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for(i in 1:ndigits){#i <- 1
		Base.b[, ndigits-i+1] <- (x %% b)
		x <- (x %/% b)
	}
	if(N ==1) Base.b[1, ] else Base.b
}
#integer.base.b(seq(0,4^3),b=4) + 1
########################################################
#Given two sequences, find all the sequences on
#the possible paths
path <- function(protein1,protein2){
  l <- length(protein1)
  U <- seq(l) # (1,2,3,4,...,l)
  diff <- protein1 - protein2
  S <- U[diff==0] #set of same sites
  D <- U[diff!=0] #set of different sites
  d <- length(D) #number of different sites
  m <- 2^(length(D))#number of sequences
  paths <- array(NA,dim=c(m,l))
  for(i in 0:(m-1)){
    pos <- binary(i,d)
    paths[i+1,] <- protein1
    paths[i+1,D[pos==1]] <- protein2[D[pos==1]]
  }
  return(paths)
}
# p1 <- c(4,3,6,8,3)
# p2 <- c(2,3,5,8,4)
# path(p1,p2)


################################################################################

#Function to choose a state given the probabilities of changing to 
#all states (probability vector)
#CONSIDER USING "SAMPLE" INSTEAD OF THIS FUNCTION
mkv <- function(odds,vec){
  vec <- vec/sum(vec)
  vec_order <- order(vec,decreasing=T)
  vec <- vec[vec_order]
  if(odds<vec[1]){
    pick <- 1
  }
  else{
    k <- 1
    while(odds > sum(vec[1:k])){
      k <- k + 1
      if(odds < sum(vec[1:k])){
        pick <- k
      }
    }
  }
  return(vec_order[pick])
}

#Given a protein, the optimal protein, and the selection vector (or constant selection)
#find the vector of rates of moving from the protein to all its neighbors
#If the length of the protein is l, then it has (m-1)*l neighbors
rate_move <- function(protein, protein_op, m, s,indep=TRUE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne)){
  l <- length(protein) #number of sites
  rates <- vector() # rates of moving to neighboring proteins, vector to return
  d1 <- pchem_d(protein,protein_op) #distance vector between given protein and optimal protein
  for(i in 1:l){ #loop over every position
    d2 <- d1 #set the distance vectors to be the same, later on, change the specific one with change
    aa <- seq(1:m)[-protein[i]] #All the m-1 aa's the site can change to
    for(j in 1:(m-1)){ #loop over the m-1 amino acids
        d2[i] <- GM[protein_op[i],aa[j]]
        rates <- c(rates,2*Ne*mu*fix(d1,d2,s,indep,formula,model,a1,a2,Phi,q,Ne))
    }
  }
  rates
}


#Simulation. Given the starting protein and the optimal protein
#t: running time of the chain
simulation <- function(protein,protein_op,t,m,s,indep=FALSE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne),record.fitness=TRUE){
  #browser()
  l <- length(protein) #number of sites
  C <- a1 + a2*l #protein production cost
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #colomn names
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  while(t_now < t){ #when current time is less than t
    rates <- rate_move(protein,protein_op,m,s,indep,formula,model,a1,a2,Phi,q,Ne,mu) #moving rates of a protein to its neighbors
    lambda <- sum(rates) # rate of staying at the same state
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait
    if(t_now < t){
      #rates <- rates/lambda #probability vector of moving
      r <- runif(1)
      index <- sample(length(rates),1,replace=T,rates) #index of the protein it moves to
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
      ftny_vec <- c(ftny_vec,Ftny_protein(path[i,seq(1:l)],protein_op,s,model)) #functionality
    }
    fitness_vec <- exp(-q*Phi*C/ftny_vec) #fitness
    path <- cbind(path,ftny_vec)
    path <- cbind(path,fitness_vec)
    colnames(path)[l+3] <- "Functionality"
    colnames(path)[l+4] <- "Fitness"
  }
  path
}

#Simulation of protein sequences of length "l" on a phylogeny "tree", 
#given the ancestral sequence "rootseq", optimal amino acid sequence "protein_op",
#selection coefficient "s", 
simTree <- function(tree, l=1000,protein_op=rep(1,l),m,s,indep=TRUE,formula=1,model=2,a1=2, a2=1, 
                    Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne),
                     bf=NULL,rootseq=NULL,ancestral=FALSE){
    if(is.null(bf)) bf = rep(1/m,m) #base frequency, randomly chosen from all states
    if(is.null(rootseq))rootseq = sample(c(1:m), l, replace=TRUE, prob=bf) #sequence at the root
    tree = ape:::reorder.phylo(tree) #oder the tree, cladwise
    edge = tree$edge #edges
    nNodes = max(edge) #number of nodes in the tree
    res = matrix(NA, nNodes, l) #result to return, amino acid sequences at all nodes
    parent <- as.integer(edge[, 1]) #parents of the edges
    child <- as.integer(edge[, 2]) #children of the edges
    root <- as.integer(parent[!match(parent, child, 0)][1])
    res[root,] <- rootseq
    tl = tree$edge.length #lengths of the edges
    #browser()
    for(i in 1:length(tl)){
        from = parent[i] 
        to = child[i]
        res[to,] = as.numeric(tail(simulation(res[from,],protein_op,tl[i],m,s,indep,
                              formula,model,a1,a2,Phi,q,Ne,mu,record.fitness=FALSE),1)[c(1:l)])
    }
    k = length(tree$tip)
    label = c(tree$tip, as.character((k+1):nNodes))
    rownames(res)=label 
    if(!ancestral)res = res[ tree$tip, , drop=FALSE]
    return(as.data.frame(res))
    #if(pt=="AA") return(phyDat.AA(as.data.frame(res), return.index=TRUE))
}
########################################################
#Find the likelihood of a tree given data
########################################################
#Likelihood of data on a tree, given the selection coefficient "s" and the optimal amino acid sequence "protein_op"
ll_indep <- function(tree,data,protein_op,m=20,s=1,rootseq=NULL,
                             formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-3,Ne=1.37e03){
  l <- length(protein_op) #number of sites in the sequence
  pr <- 0
  for(i in 1:l){
    llh = ll_site(tree,data[,i],protein_op[i],m,s,rootseq,formula,model,a1,a2,Phi,q,Ne)
    pr = pr - log(llh)
  }
  pr
}

ll_site <- function(tree,data,optimal,m=20,s=1,rootseq=NULL,
                            formula=1,model=2,a1=2,a2=1,Phi=0.5,q=4e-3,Ne=1.37e03){
    
    if(is.null(rootseq)) rootseq = rep(1/m,m) #base frequency, randomly chosen from all states
    Q = mat_gen_indep(optimal,m,s,formula,model,a1,a2,Phi,q,Ne) #transition rate matrix for the site, given the optimal aa
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
        for(j in 1:m){
          #find the likelihood of a node when likelihood of both its descendents is found --- pruning
          probvec[from,j] = (expm(Q*t_left)[j,] %*% probvec[to[1],])*(expm(Q*t_right)[j,] %*% probvec[to[2],]) 
        }
    }
    browser()
    probvec[root,] %*% rootseq
}

