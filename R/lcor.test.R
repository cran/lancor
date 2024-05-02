lcor.test = function(x, y = NULL, type = c("rank", "linear"), nperm = 999,
                     method = c("permutation", "asymptotic", "symmetric")) {
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (is.matrix(x)) {
    y = x[,2]
    x = x[,1]
  }
  type = match.arg(type)
  if (!(type %in% c("rank", "linear")))
    stop("type for lcor can only be rank or linear")
  #
  if (type == "rank") {
    lc = lcor(x, y, type="rank")
    n = length(x)
    method = match.arg(method)
    if (!(method %in% c("asymptotic", "permutation")))
      stop("method for rank test can only be permutation or asymptotic")
    if (method == "asymptotic") { #asymptotic theory
      ts = sqrt(n) * lc
      pval = 1 - (2*pnorm(ts) - 1)^2
      return(list(lcor = lc, pval = pval))
    }
    if (method == "permutation") {
      ts = lc
      if (n <= 6) {
        nperm = arrangements::npermutations(n) #use all permutations
        perm = arrangements::permutations(x)
      } else {
        perm = arrangements::permutations(x, nsample = nperm)
      }
      tp = rep(0, nperm)
      for (i in 1:nperm) {
        tp[i] = lcor( perm[i,], y, type="rank")
      }
      pval = (sum(tp > ts) + 1) / (nperm + 1)
      return(list(lcor = lc, pval = pval))
    }
  }
  #
  if (type == "linear") {
    lc = lcor(x, y, type="linear")
    n = length(x)
    if (!(method %in% c("permutation", "asymptotic", "symmetric")))
      stop("method for test can only be permutation, asymptotic or symmetric")
    if (method == "symmetric") { #asymptotic theory with vanishing third moments
      ts = sqrt(n) * lc
      pval = 1 - (2*pnorm(ts) - 1)^2
      return(list(lcor = lc, pval = pval))
    }
    if (method == "asymptotic") { #asymptotic theory with arbitrary third moments
      ts = sqrt(n) * lc
      x = sqrt(n/(n-1)) * scale(x) # factor ensures that mean(x^2)=1, hence e40>1
      y = sqrt(n/(n-1)) * scale(y)
      e30 = mean(x^3)
      e40 = mean(x^4)
      e03 = mean(y^3)
      e04 = mean(y^4)
      tau = e30*e03 / sqrt( (e40-1) * (e04-1) )
      alpha1 = sqrt( (1-tau) / (1+tau) )
      alpha2 = sqrt( (1+tau) / (1-tau) )
      pval = 1 - 2 * ( sn::psn(ts,0,1,alpha1) - sn::psn(0,0,1,alpha1)
                       + sn::psn(-ts,0,1,alpha2) - sn::psn(0,0,1,alpha2) )
      return(list(lcor = lc, pval = pval))
    }
    if (method == "permutation") {
      ts = lc
      if (n <= 6) {
        nperm = arrangements::npermutations(n) #use all permutations
        perm = arrangements::permutations(x)
      } else {
        perm = arrangements::permutations(x, nsample = nperm)
      }
      tp = rep(0, nperm)
      for (i in 1:nperm) {
        tp[i] = lcor( perm[i,], y, type="linear")
      }
      pval = (sum(tp > ts) + 1) / (nperm + 1)
      return(list(lcor = lc, pval = pval))
    }
  }
}
