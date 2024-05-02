ace.test = function(x, y = NULL, nperm = 999) { #permutation test with ace
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (is.matrix(x)) {
    y = x[,2]
    x = x[,1]
  }
  a = acepack::ace(x, y)
  ace.cor =  as.vector( cor(a$tx, a$ty) )
  n = length(x)

  ts = ace.cor
  if (n <= 6) {
    nperm = arrangements::npermutations(n) #use all permutations
    perm = arrangements::permutations(x)
  } else {
    perm = arrangements::permutations(x, nsample = nperm)
  }
  tp = rep(0, nperm)
  for (i in 1:nperm) {
    a = acepack::ace(perm[i,], y)
    tp[i] = as.vector( cor(a$tx, a$ty) )
  }
  pval = (sum(tp > ts) + 1) / (nperm + 1)
  return(list(ace = ace.cor, pval = pval))
}
