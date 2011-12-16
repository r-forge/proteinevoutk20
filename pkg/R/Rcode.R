##################################
#physiochemical distance matrix#
##################################
#notice that the proteins in matrix rows and columns are not in the same order
library("Matrix")
#parameter values
a1 <- 2
a2 <- 1
Phi <- 0.5
q <- 4*10^(-7)
Ne <- 1.37*10^7
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
GM <- GM/100 #make the mean 1 instead of 100

#find the physiochemical distance vector between two proteins
pchem_d <- function(protein1, protein2){
  l <- length(protein1)
  d <- vector(length=l)
  for(i in 1:l){
    d[i] <- GM[protein1[i],protein2[i]]
  }
  d
}
#############################################################################
#Functions to find the means for categories of discrete gamma distribution
#x*dgamma(x)
xgamma <- function(x,shape,scale){
  return(x*dgamma(x,shape,rate=scale))
}

#Given the shape and scale of the gamma distribution,
#the number of categories of discretization.
#Find the mean of portion of the gamma distribution falling in 
#each category
#THIS FUNCTION WILL RETURN ERROR WHEN THE SHAPE PARAMETER IS TOO SMALL
#BECAUSE THE INTEGRATION VALUE WILL BE NON-FINITE!!!!
dis_gamma <- function(shape,scale,category=4){
  #browser()
  z <- vector(length=(category+1))
  x <- vector(length=category)
  z[1] <- 0
  z[category+1] <- Inf
  for(i in 2:category){
    z[i] <- qgamma((i-1)/category,shape=shape,rate=scale)
  }
  for(j in 1:category){
    x[j] <- integrate(xgamma,shape=shape,scale=scale,lower=z[j],upper=z[j+1])$value/(integrate(dgamma,shape=shape,rate=scale,lower=z[j],upper=z[j+1])$value)
  }
  return(x)
}

#means <- dis_gamma(sh,sc)
#################################################################
#Functionality function, arguments: 
#d: distance between the optimal and observed amino acids, a vector of length n
#s: seleciton strength coeffecient vector, comes from a probability distribution, e.g. exponential or gamma
#model: model number. 1: exponential 2: parabolic
Ftny <- function(d, s, model=1){  #default: n=15, exponential distribution with rate 1
  #browser()
  if(model==1){ #model 1
    result <- prod(exp(-d*s)) #the function
  }
  else if(model==2){ #model 2
    result <- prod(1/(1+d*s)) 
  }
  return(result)
} #d and s are given

#fitness function
#q, Phi: constants
#C: cost function, linearly function of the length of the protein
#C_n = a1 + a2*n, cost --  linear function of the length

#fixation probability function from protein 1 to protein 2,
#which is also the instantaneous rate from protein1 to protein2
#parameters: d1, d2: distances from the optimal protein, vectors
# s: selection coefficients for all sites
# a1, a2: parameters in cost function C_n
#Ne: effective population size
#formula number: 1: Sela & Hirsh's formula (default); 2: classic formula
fix <- function(d1, d2, s, formula=1, model=1){
  #browser()
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
      ftny_S <- Ftny(d1[S],s[S],model) #same terms in the product for both proteins
      ftny_D1 <- Ftny(d1[D],s[D],model) #different terms in the product
      ftny_D2 <- Ftny(d2[D],s[D],model)
      if(ftny_D1==ftny_D2){
        return(1/(2*Ne))
      }
      else{
        fit_ratio <- exp(-q*Phi*C_n*(1/ftny_S)*((1/ftny_D1)-(1/ftny_D2)))
      }
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
        result <- (1-fit_ratio)/(1-fit_ratio^(2*Ne))
      }
      else if(formula==2){
        se <- 1/fit_ratio - 1 # s=(w1-w2)/w2
        #result <- (1-exp(-2*se))/(1-exp(-4*Ne*se))
        result <- (1-exp(-se))/(1-exp(-2*Ne*se))
      }
      else{
        print("Error! Please choose 1 or 2 as formula number. 1: Hirsh's; 2: classic pop gen's")
      } 
  }
  return(result)
}

#Assign values to different parameters

######################################################

#Function to generate the transition matrix (20*20 for now) 
#given the optimal amino acid. This is only for one site (amino acid)
#parameter: aa_op: optimal amino acid at this site
#s: selection coefficient
mat_gen_indep <- function(aa_op,s,formula=1,model=1){
  mat <- matrix(0,nrow=20,ncol=20) #set diagonal entries to be 0 at first
  for(i in 1:19){
    for(j in (i+1):20){
      d1 <- GM[i,aa_op] #distance between aa_i and aa_op
      d2 <- GM[j,aa_op]
      mat[i,j] <- fix(d1,d2,s,formula,model) #fixation prob -> transition rate
      mat[j,i] <- fix(d2,d1,s,formula,model) #symmetric entry
    }#end for j
    mat[i,i] <- -(sum(mat[i,])) #diagonal entry; sum of each row is 0
  }#end for i
  mat[20,20] <- -sum(mat[20,])
  return(mat)
}

#arr: the protein sequences; each row represents a protein considered in the matrix
#op: optimal amino acid sequences
#s: selection coefficients
#formula, model
#########################################################################3
mat_gen_dep <- function(arr,protein_op,s,formula=1,model=1){
  #browser()
  m <- dim(arr)[1] #number of proteins to consider
  n <- dim(arr)[2] #number of sites in each protein
  d1 <- rep(0,n)
  d2 <- d1
  mat <- matrix(0,nrow=m,ncol=m) #matrix to return, every entry is initiated to be 0
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      d1 <- pchem(arr[i,],protein_op)
      d1 <- pchem(arr[j,],protein_op)
      mat[i,j] <- fix(d1,d2,s,formula,model) #fixation prob -> transition rate
      mat[j,i] <- fix(d2,d1,s,formula,model) #symmetric entry
    }
    mat[i,i] <- -sum(mat[i,])
  }
  mat[m,m] <- -sum(mat[m,])
  return(mat)
}
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
####################################################################################################
#Given two protein sequences of the same lengths, calculated the probability (-loglikelihood) of 
#going from protein1 to protein2 using site independent method.
#Calculate the probability of going from one amino acid to another, then
#multiply them together to get the total probability.
#Parameter: branch length t. There is only one branch and two nodes in this simplest case
prob_indep_gamma <- function(protein1, protein2, protein_op, t,shape,category=4,formula=1,model=1){
  #browser()
  l <- length(protein1) #length of protein
  means <- dis_gamma(shape,shape,category) #means for discrete gamma distribution
  pr <- 0
  for(i in 1:l){ # for each amino acid in the protein, calculate the probability of transition
    p_t <- vector() #store the probabilities for all categories
    for(j in 1:category){
        ma <- mat_gen_indep(protein_op[i],means[j],formula, model) #transition matrix Q
        Pt <- expm(ma*t) #P(t) = exp(Qt)
        p_t <- c(p_t,Pt[protein1[i],protein2[i]])#the entry corresponding to the amino acids in the sequences
    }
    p <- mean(p_t) #f(x)=(1/k)*sum(f(x|means))
    pr <- pr-log(p)
  }  
  return(pr)
}

prob_dep_gamma <- function(protein1,protein2,protein_op,t,shape,category=4,formula=1,model=1){
  Arr <- path(protein1,protein2) # all the sequences lying on the path from protein 1 to protein 2
  m <- dim(Arr)[1]
  l <- length(protein1) #number of amino acids
  means <- dis_gamma(shape,shape,category) #means in discrete gamma distribution for each category
  means_index <- integer.base.b(seq(0,category^l-1),category)+1 #all possible combinations of (discrete) gamma means
  index <- dim(means_index)[1] #number of assignments of means to each site (s)
  pr <- vector(length=index)
  for(i in 1:index){
    s <- means[means_index[i,]]
    ma <- mat_gen_dep(Arr,protein_op,s,formula,model)
    Pt <- expm(ma*t)
    pr <- c(pr,Pt[1,m]) #first row last column entry of P(t)
  }
  p <- mean(pr)
  return(-log(p))
}


prob_gamma <- function(protein1, protein2, protein_op, t, shape, indep=TRUE,category=4,formula=1,model=1){
  if(indep){ #site independent case
    prob_indep_gamma(protein1,protein2,protein_op,t,shape,category,formula,model)
  }
  else{ #site dependent case
    prob_dep_gamma(protein1,protein2,protein_op,t,shape,category,formula,model)
  }
  return(pr)
}
###################################################################################################
#Fixed rate for s -> likelihood
#Independent-site version
prob_indep_fix <- function(protein1,protein2,protein_op,t,s=1,formula=1,model=1){
  l <- length(protein1) #length of protein
  pr <- 0
  for(i in 1:l){ # for each amino acid in the protein, calculate the probability of transition
    ma <- mat_gen_indep(protein_op[i],s,formula, model) #transition matrix Q
    Pt <- expm(ma*t) #P(t) = exp(Qt)
    p <- Pt[protein1[i],protein2[i]]#the entry corresponding to the amino acids in the sequences
    pr <- pr-log(p)
  }  
  return(pr)
}

#Dependent-site version
prob_dep_fix <- function(protein1,protein2,protein_op,t,s=1,formula=1,model=1){
  Arr <- path(protein1,protein2) # all the sequences lying on the path from protein 1 to protein 2
  m <- dim(Arr)[1] # number of amino acid sequences on the parsimony path
  l <- length(protein1) #number of amino acids in the protein
  sv <- rep(s,l) #s is the same for all sites
  ma <- mat_gen_dep(Arr,protein_op,sv,formula,model) #transition rate matrix for the path Q
  Pt <- expm(ma*t) # P(t) = Exp(Q*t)
  p <- Pt[1,m] #first row last column entry of P(t)
  return(-log(p))
}

prob_fix <- function(para, protein1, protein2, protein_op, indep=TRUE,formula=1,model=1){
  t <- para[1]
  s <- para[2]
  if(indep){ #site independent case
    prob_indep_fix(protein1,protein2,protein_op,t,s,formula,model)
  }
  else{ #site dependent case
    prob_dep_fix(protein1,protein2,protein_op,t,s,formula,model)
  }
}



#Function to choose a state given the probabilities of changing to 
#all states (probability vector)
mkv <- function(odds,vec){
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
  return(pick)
}

rate_move_dep <- function(protein, protein_op, s,formula=1,model=1){
  l <- length(protein)
  sv <- rep(s,l)
  rates <- vector() #store the rates of moving to other states
  d1 <- vector(length=l) #distance from optimal aa's
  for(i in 1:l){ 
    d1[i] <- GM[protein[i],protein_op[i]]
  }
  #browser()
  for(i in 1:l){
    d2 <- d1
    aa <- seq(1:20)[-protein[i]]
    for(j in 1:19){
      d2[i] <- GM[protein_op[i],aa[j]]
      rates <- c(rates,fix(d1,d2,sv,formula,model))
    }
  }
  rates
}
#rate_move_dep(c(1,2),c(3,3),c(1,1))

simulation_dep <- function(protein,protein_op,t,s=1,formula=1,model=1,record.fitness=TRUE){
  #browser()
  l <- length(protein) #number of sites
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #colomn names
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  while(t_now < t){ #when current time is less than t
    rates <- rate_move_dep(protein,protein_op,s,formula,model) #moving rates of a protein to its neighbors
    lambda <- sum(rates) 
    t_wait <- rexp(1,lambda) #waiting time
    t_now <- t_now + t_wait
    if(t_now < t){
      rates <- rates/lambda #probability vector of moving
      index <- mkv(runif(1),rates) #index of the protein it moves to
      pos <- (index-1) %/% 19 + 1 #index of the site that changes
      rmd <- (index-1) %% 19 + 1 # the amino acid that site changes to
      protein[pos] <- aa[-protein[pos]][rmd]
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  if(record.fitness){
    fitness_vec <- NULL
    for(i in 1:dim(path)[1]){
      fitness_vec <- c(fitness_vec,Ftny(pchem_d(path[i,seq(1:l)],protein_op),model))
    }
    path <- cbind(path,fitness_vec)
  }
  path
}
#simulation_dep(rep(3,3),seq(1:3),10^7,1)

#simulation for site-independent case
#Multiple sites(amino acids) in the protein
simulation_indep <- function(protein, protein_op, t, s=1,formula=1,model=1,record.fitness=TRUE){
  l <- length(protein) #number of sites(amino acids)
  t_now <- 0 #time until now
  arr <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(arr) <- colnames(arr,do.NULL = FALSE, prefix = "Site.")
  colnames(arr)[l+1] <- "Time Now"
  colnames(arr)[l+2] <- "Waiting Time"
  while(t_now < t){ #while the time does not exceed the simulation time
    t_sim <- NULL # simulated waiting times
    for(i in 1:l){
      rate <- sum(find_rate(protein[i],protein_op[i],s,formula, model)) # rates of exponential distributions for waiting times
      t_sim <- c(t_sim,rexp(1,rate)) #waiting times
    }
    t_wait <- min(t_sim) #shortest waiting time
    pos <- which.min(t_sim) #site where the shortest waiting time occurs
    t_now <- t_now + t_wait #time after waiting (up to now)
    if(t_now < t){
      #probability vector of moving to other aa's at "pos"
      move_prob <- find_rate(protein[pos],protein_op[pos],s,formula,model)/sum(find_rate(protein[pos],protein_op[pos],s,formula,model))
      protein[pos] <- mkv(runif(1),move_prob)#the amino acid it moves to
      arr <- rbind(arr,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      arr <- rbind(arr,c(protein,t,NA)) #record this moving step
    }
  }
    if(record.fitness){
    fitness_vec <- NULL
    for(i in 1:dim(arr)[1]){
      fitness_vec <- c(fitness_vec,Ftny(pchem_d(arr[i,seq(1:l)],protein_op),model))
    }
    arr <- cbind(arr,fitness_vec)
  }
  arr
}
#simulation_indep(rep(2,3),seq(4:6),10^6)
#find the instantaneous rates of moving from "protein" to another 
# amino acid, with optimal amino acid "protein_op"
# Note: the rate of not moving is set to 0, here, the sum of the vector is not equal to 0
find_rate <- function(protein, protein_op, s=1, formula=1, model=1){
  aa <- seq(1:20)
  rate <- vector()
  for(i in 1:length(aa)){
    if(i==protein){
      rate <- c(rate,0)
    }
    else{
    rate <- c(rate,fix(GM[protein,protein_op],GM[aa[i],protein_op],s,formula,model))
    }
  }
  rate
}
#optim(c(10,0.1),prob_fix,protein1=p1,protein2=p2,protein_op=protein_op)
#optim(c(10,1),prob_fix,protein1=p1,protein2=p2,protein_op=protein_op)
################################################################
#Optimization: maximum likelihood estimates for parameters
################################################################

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
####################################################################################
Ftny_sim <- function(d, s, model=1, dis=0, para=1){
  #browser()
  #different distributions where s comes from
  n <- length(d) #length of the protein as an amino acid sequence
  if(dis==0){
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
} #This is used to do simulation, not the calculation