##################################
#physiochemical distance matrix#
##################################
#notice that the proteins in matrix rows and columns are not in the same order
GranthamMatrix <- read.csv("GranthamMatrix.csv",
                           header=TRUE, sep=",",row.names=1)
#View(GranthamMatrix) #check the matrix
GM <- GranthamMatrix #store it in a different matrix and manipulate it
#complete the upper triangle matrix to whole, with diagonals equal to 0
GM[1,1] <- 0
for(i in 2:20){
  GM[i,i] <- 0
  for(j in 1:(i-1)){
    GM[i,j] <- GM[j,i]
  }
}
GM <- GM/100 #make the mean 1 instead of 100
############################################
#Functionality function, arguments: 
#d: distance between the optimal and observed amino acids, a vector of length n
#s: seleciton strength coeffecient vector, comes from a probability distribution, e.g. exponential or gamma
#n: length of the protein
#dis: distribution function from which s is drawn. 1: exponential, 2: gamma, 3: uniform(0,1)
#para: parameter for the distribution. Rate for exponential, shape for gamma
#model: model number. 1: exponential 2: parabolic
Ftny <- function(d, s, model=1, dis=0, para=1){  #default: n=15, exponential distribution with rate 1
  #browser()
  #different distributions where s comes from
  n <- length(d) #length of the protein as an amino acid sequence
  if(dis==0){
    s <- rep(s,n)
  }
  else if(dis==1)
    {s <- rexp(n,para)} #random generation of n numbers from exponential distribution
  else if(dis==2) #gamma distribution
    {s <- rgamma(n,para)}
  else if(dis==3) # uniform distribution
    {s <- runif(n)}

  #different models for functionality function
  if(model==1){ #model 1
    result <- prod(exp(-d*s)) #the function
  }
  else if(model==2){ #model 2
    result <- prod(1/(1+d*s)) 
  }
  return(result)
}

#fitness function
#q, Phi: constants
#C: cost function, linearly function of the length of the protein
#ftn: functionality value
#C_n = a1 + a2*n, cost --  linear function of the length
#n: length of the protein

#fixation probability function from protein 1 to protein 2,
#which is also the instantaneous rate from protein1 to protein2
#parameters: fit1, fit2: fitnesses of protein 1 and 2
# n: length of protein
# a1, a2: parameters in cost function C_n
#Ne: effective population size
#formula number: 1: Sela & Hirsh's formula (default); 2: classic formula

#####
#try to put it in a further calculated form, with same terms for ftny's pulle out
#as common factors. To do this, have to recognize which ONE is the one that's different
#in the two proteins
#####
fix <- function(d1, d2, a1, a2, Ne, formula=1,s=0.8,model=1,dis=0,para=1){
  #browser()
  n <- length(d1)
  C_n <- a1 + a2*n
  if(n==1){ #if there is only 1 amino acid in the sequence
    ftny1 <- Ftny(d1,s,model,dis,para)
    ftny2 <- Ftny(d2,s,model,dis,para)
    if(ftny1==ftny2){ #When the fitnesses are the same, neutral case, pure drift
      return(1/(2*Ne))
    }
    else{
      fit_ratio <- exp(-q*Phi*C_n*((1/ftny2)-(1/ftny1))) #w2/w1
    }
  }
  else{ # there are more than 1 amino acids in the protein
    index <- seq(1,n)
    diff <- d1 - d2
    S <- index[diff==0]
    D <- index[diff!=0]
    if(length(D) > 1){ #if there are more than 1 site with different amino acids
      return(0)
    }
    else if(length(D)==1){ #if there is only one site that is different
      ftny_S <- Ftny(d1[S],s,model,dis,para)
      ftny_D1 <- Ftny(d1[D],s,model,dis,para)
      ftny_D2 <- Ftny(d2[D],s,model,dis,para)
      if(ftny_D1==ftny_D2){
        return(1/(2*Ne))
      }
      else{
        fit_ratio <- exp(-q*Phi*C_n*(1/ftny_S)*((1/ftny_D2)-(1/ftny_D1)))
      }
    }
  }
  #######
  #different fixation prob formula
  #######
      if(formula==1){
        result <- (1-fit_ratio)/(1-fit_ratio^(2*Ne))
      }
      else if(formula==2){
        se <- fit_ratio - 1 # s=(w2-w1)/w1
        result <- (1-exp(-2*se))/(1-exp(-4*Ne*se))
      }
      else{
        print("Error! Please choose 1 or 2 as formula number. 1: Hirsh's; 2: classic pop gen's")
      } 

  return(result)
}

#Assign values to different parameters

######################################################

#Function to generate the transition matrix (20*20 for now) 
#given the optimal amino acid
#parameter: op(aa_op): optimal amino acid at this site
mat_gen_indep <- function(op){
  mat <- matrix(0,nrow=20,ncol=20) #set diagonal entries to be 0 at first
  for(i in 1:19){
    for(j in (i+1):20){
      d1 <- GM[i,op] #distance between aa_i and aa_op
      d2 <- GM[j,op]
      mat[i,j] <- fix(d1,d2,a1,a2,Ne,formula=1) #fixation prob -> transition rate
      mat[j,i] <- fix(d2,d1,a1,a2,Ne,formula=1) #symmetric entry
    }#end for j
    mat[i,i] <- -(sum(mat[i,])) #diagonal entry; sum of each row is 0
  }#end for i
  mat[20,20] <- -sum(mat[20,])
  return(mat)
}

#take care of the NaN's in the matrix
#########################################################################3
matgen_dep <- function(arr,op){
  m <- dim(arr)[1]
  mat <- matrix(0,dim=c(m,m))
  for(i in 1:(m-1)){
    mat[i,i] <- 0
    for(j in (i+1):m){
      d1 <- rep(0,n)
      d2 <- d1
      for(k in 1:n)
        d1[k] <- GM[arr[i,k],op]#distance vectors
        d2[k] <- GM[arr[j,k],op]
        mat[i,j] <- fix(d1,d2,a1,a2,Ne,formula=1) #fixation prob -> transition rate
        mat[i,j] <- fix(d2,d1,a1,a2,Ne,formula=1) #symmetric entry
    }
    mat[i,i] <- -sum(mat[i,])
  }
  mat[m,m] <- -sum(mat[m,])
  return(mat)
}
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


#Given two sequences, find all the sequences on
#the possible paths
path <- function(p1,p2){
  l <- length(p1)
  U <- seq(l) # (1,2,3,4,...,l)
  diff <- p1 - p2
  S <- U[diff==0] #set of same sites
  D <- U[diff!=0] #set of different sites
  d <- length(D) #number of different sites
  m <- 2^(length(D))#number of sequences
  paths <- array(NA,dim=c(m,l))
  for(i in 0:(m-1)){
    pos <- binary(i,d)
    paths[i+1,] <- p1
    paths[i+1,D[pos==1]] <- p2[D[pos==1]]
  }
  return(paths)
}
p1 <- c(4,3,6,8,3)
p2 <- c(2,3,5,8,4)
path(p1,p2)

#Given two protein sequences of the same lengths, calculated the probability (likelihood) of 
#going from protein1 to protein2 using site independent method.
#Calculate the probability of going from one amino acid to another, then
#multiply them together to get the total probability.
#Parameter: branch length t. There is only one branch and two nodes in this simplest case
prob <- function(prn1, prn2, prn_op, t, indep=TRUE){
  l <- length(prn1) #length of protein
  if(indep){ #site independent case
    pr <- 1
    for(i in 1:l){
      ma <- mat_gen_indep(prn_op[i]) #transition matrix Q
      Pt <- expm(ma*t) #P(t) = exp(Qt)
      pr <- pr*pt[prn1[i],prn2[2]] #the entry corresponding to the amino acids in the sequences
    }
  }
  else{ #site dependent case
    paths <- path(prn1,prn2)
    ma <- mat_gen_dep(paths,prn_op)
    Pt <- expm(ma*t)
    pr <- Pt[1,dim(Pt)[2]] #first row last column entry of P(t)
  }
  return(pr)
}

################################################################
#Optimization: maximum likelihood estimates for parameters
################################################################
  #parameter values
Phi <- 0.25
q <- 4*10^(-7)
Ne <- 1.37*10^7
########################################################################
########################################################################
#   pr1 <- sample(1:20,n,replace=T) #amino acids in protein1
#   pr2 <- pr1 #amino acids in protein2
#   pos <- sample(1:n,1) # position where pr2 has different aa from pr1
#   pr2[pos] <- sample(1:20,1) #assign a different amino acid at this position
#   pr0 <- sample(1:20,n,replace=T) #amino acids in optimal protein
#   d1 <- vector(length=n) #distance vector from pr1 to pr0
#   d2 <- d1 #distance vector from pr2 to pr0
#   for(i in 1:n){
#     d1[i] <- GM[pr1[i],pr0[i]]
#     d2[i] <- GM[pr2[i],pr0[i]]
#   }

