gos <- lapply(1:length(GGList), function(x) GGList[[x]]$order)
sp <- scatterplot3d(sw_max[1,gos[[5]]],sw_max[2,gos[[5]]],log(sw_max[3,gos[[5]]]),
                    xlab="g", ylab="beta",zlab="log(gamma)")
sp$points3d(sw_max[1,gos[[2]]],sw_max[2,gos[[2]]],log(sw_max[3,gos[[2]]]),col="blue")
sp$points3d(sw_max[1,gos[[3]]],sw_max[2,gos[[3]]],log(sw_max[3,gos[[3]]]),col="red")
sp$points3d(sw_max[1,gos[[4]]],sw_max[2,gos[[4]]],log(sw_max[3,gos[[4]]]),col="green")

sp1 <- scatterplot3d(sw_maj[1,gos[[1]]],sw_maj[2,gos[[1]]],log(sw_maj[3,gos[[1]]]),
                     xlab="g", ylab="beta",zlab="log(gamma)")
sp1$points3d(sw_maj[1,gos[[2]]],sw_maj[2,gos[[2]]],log(sw_maj[3,gos[[2]]]),col="blue")
sp1$points3d(sw_maj[1,gos[[3]]],sw_maj[2,gos[[3]]],log(sw_maj[3,gos[[3]]]),col="red")
sp1$points3d(sw_maj[1,gos[[4]]],sw_maj[2,gos[[4]]],log(sw_maj[3,gos[[4]]]),col="green")


var_s_max <- sapply(1:length(gos), function(x) var(s_max[gos[[x]]]))
var_s_maj <- sapply(1:length(gos), function(x) var(s_maj[gos[[x]]]))

#for each gene, find the amino acids that have nonzero weights of being optimal
# using opw_opw1 (opw_opw has one gene 25 that doesn't have a result)
opAA <- vector(mode="list",length=106)
for(i in 1:106){
  index <- which(opw_opw1[,i]!=0)
  value <- opw_opw1[index,i]
  opAA[[i]]$AA <- index[order(value,decreasing=TRUE)]
  opAA[[i]]$weights <- value[order(value,decreasing=TRUE)]
}

s_list <- vector(mode="list",length=length(gos))
for(i in 1:length(gos)){
  s_list[[i]] <- s_max[gos[[i]]]
}
beta_list <- vector(mode="list",length=length(gos))
for(i in 1:length(gos)){
  beta_list[[i]] <- GM_max[2,gos[[i]]]
}
gamma_list <- vector(mode="list",length=length(gos))
for(i in 1:length(gos)){
  gamma_list[[i]] <- GM_max[3,gos[[i]]]
}