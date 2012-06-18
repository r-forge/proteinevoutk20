> simSeq1 <- 
function (tree, l = 1000, Q = NULL, bf = NULL, rootseq = NULL, 
    type = "DNA", model = "USER", levels = NULL, rate = 1, ancestral = FALSE) 
{
	browser()
    pt <- match.arg(type, c("DNA", "AA", "USER"))
    if (pt == "DNA") 
        levels <- c("a", "c", "g", "t")
    if (pt == "AA") 
        levels <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", 
            "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", 
            "v")
    if (pt == "USER") 
        if (is.null(levels)) 
            stop("levels have to be supplied if type is USER")
    lbf = length(levels) #number of states for the characters
    if (type == "AA" & !is.null(model)) {
        model <- match.arg(model, c("USER", "WAG", "JTT", "LG", 
            "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        if (model != "USER") 
            getModelAA(model, bf = is.null(bf), Q = is.null(Q))
    }
    if (is.null(bf)) 
        bf = rep(1/lbf, lbf) #base frequency?
    if (is.null(Q)) 
        Q = rep(1, lbf * (lbf - 1)/2)
    if (is.matrix(Q)) 
        Q = Q[lower.tri(Q)]
    eig = edQt(Q, bf) #eigenvalues, Q, and Q inverse
    m = length(levels)
    if (is.null(rootseq)) 
        rootseq = sample(levels, l, replace = TRUE, prob = bf)
    tree = ape:::reorder.phylo(tree)
    edge = tree$edge
    nNodes = max(edge)
    res = matrix(NA, l, nNodes)
    parent <- as.integer(edge[, 1])
    child <- as.integer(edge[, 2])
    root <- as.integer(parent[!match(parent, child, 0)][1])
    res[, root] = rootseq
    tl = tree$edge.length
    for (i in 1:length(tl)) {
        from = parent[i]
        to = child[i]
        P = getP(tl[i], eig, rate)[[1]]
        for (j in 1:m) {
            ind = res[, from] == levels[j]
            res[ind, to] = sample(levels, sum(ind), replace = TRUE, 
                prob = P[, j])
        }
    }
    k = length(tree$tip)
    label = c(tree$tip, as.character((k + 1):nNodes))
    colnames(res) = label
    if (!ancestral) 
        res = res[, tree$tip, drop = FALSE]
    if (pt == "DNA") 
        return(phyDat.DNA(as.data.frame(res), return.index = TRUE))
    if (pt == "AA") 
        return(phyDat.AA(as.data.frame(res), return.index = TRUE))
    if (pt == "USER") 
        return(phyDat.default(as.data.frame(res), levels = levels, 
            return.index = TRUE))
}
<environment: namespace:phangorn>