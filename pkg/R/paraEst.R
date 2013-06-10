# nchar=c(20, 200, 2000)
# true.parameters<-somehow
# for (nchar.index in sequence(length(nchar))) {
#   optimal.aa<-rep(c(1:20), nchar/20)
#   starting.aa<-optimal.aa
#   simulated.root<-single.sim(starting.aa, time=really long)
#   sims.from.tree <- single.sim(simulated.root, phy)
#   parameter.estimates <- estimate(sims.from.tree, phy)
#   plot(parameter.estimates, true.parameters)
# }
nchars = c(20,100,200,400,600,800,1000,2000)
ind = Sys.getenv("SGE_TASK_ID")
ind <- as.integer(ind)
load("~/BackupProEvo/Newton/rokas/rootMax/gene1.RData")
nchar = nchars[ind]
opaa <- rep(1:20,times=nchar/20)
nchar <- length(opaa)
start.aa <- opaa
sim.long <- simulation1(protein=opaa,protein_op=opaa,t=30,matall=res_op$Qall)
## check and see if the equilibrium has been reached
sim_info <- sim.info(sim.long$path,opaa=opaa,s=s,beta=beta,gamma=gamma)
plot(sim_info$fty~sim_info$t)

s <- res_op$s
beta <- res_op$GMweights[2]
gamma <- res_op$GMweights[3]
Q <- res_op$Q
tree <- res_op$tree
sim.root <- as.integer(tail(sim.long$path,1)[1:nchar])
sim <- simTree(tree=res_op$tree,protein_op=opaa,matall=res_op$Qall,rootseq=sim.root,scale.vec=rep(1,20))
dataphy <- phyDat(sim$data,type="AA")
res <- mllm1(data=dataphy,tree=res_op$tree,s=s,beta=be,gamma=ga,ancestral="max")
res$ll$loglik
## estimate parameters from simulated data
res_op_sim <- optim.mllm1(res,optQ=TRUE,optBranch=TRUE,optsWeight=TRUE,optscale=FALSE,optOpw=FALSE,
                      control=list(epsilon=1e-08,hmaxit=2,htrace=TRUE,print_level=0,maxeval="50"),ancestral="max")
res_op_sim$ll$loglik
