tree <- unroot(ROKAS_TREE)
if (is.null(attr(tree, "order")) || attr(tree, "order") ==
      "cladewise")
  tree <- reorder(tree, "postorder")
nr <- as.integer(attr(data, "nr")) #number of different site patterns                                                                                                        
nc <- as.integer(attr(data, "nc")) #number of states                                                                                                                         
nTips <- as.integer(length(tree$tip.label))
node <- tree$edge[, 1]
edge <- tree$edge[, 2]
m = length(edge) + 1

el <- as.double(tree$edge.length)

#    node = as.integer(node - nTips - 1L)                                                                                                                                        
node = as.integer(node - min(node))
edge = as.integer(edge - 1L)

contrast = attr(data, "contrast")
nco = as.integer(dim(contrast)[1])

root <- as.integer(getRoot(tree))
