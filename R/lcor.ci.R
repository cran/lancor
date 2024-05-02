lcor.ci = function(x, y = NULL, conf.level = 0.95, type = c("rank", "linear"),
                   con = TRUE, R = 1000, method = c("plugin", "boot", "pretest")) {
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!is.matrix(x)) x = cbind(x, y)
  xx = x
  type = match.arg(type)
  if (!(type %in% c("rank", "linear")))
    stop("type for lcor can only be rank or linear")
  requireNamespace("boot", quietly = TRUE)
  #
  if (type == "rank") {
    lc.res = lcor.comp(xx, type="rank")
    rho1 = lc.res[1]
    rho2 = lc.res[2]
    lc = lc.res[3]
    n = dim(xx)[1]
    alpha = 1 - conf.level
    #compute bootstrap covariance matrix of (rho1,rho2)
    Sigma =  n * cov( boot::boot(xx, statistic = function(xx, id) lcor.comp( xx[id, ], type="rank"), R = R)$t[,1:2] )
    #case |rho1| != |rho2|
    if (abs(rho1) >= abs(rho2)) {
      se = sqrt(Sigma[1,1])
    } else {
      se = sqrt(Sigma[2,2])
    }
    z = qnorm(1 - alpha/2)
    l1 = max( lc - z * se / sqrt(n), 0)
    u1 = min( lc + z * se / sqrt(n), 1)
    #case |rho1|=|rho2|>0
    F = function(x, Sigma, p) {
      sigma1 = sqrt(Sigma[1,1])
      sigma2 = sqrt(Sigma[2,2])
      tau = Sigma[1,2] / ( sigma1 * sigma2 )
      alpha1 = (sigma1/sigma2-tau) / sqrt(1-tau^2)
      alpha2 = (sigma2/sigma1-tau) / sqrt(1-tau^2)
      0.5 * ( sn::psn(x,0,sigma1,alpha1) + sn::psn(x,0,sigma2,alpha2) ) - p
    }
    si = 5 * max(Sigma[1,1],Sigma[2,2]) #guess for search interval
    c.u = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = alpha/2)$root
    c.o = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = 1-alpha/2)$root
    l2 = max( lc - c.o / sqrt(n), 0)
    u2 = min( lc - c.u / sqrt(n), 1)
    if (con==TRUE) return(setNames(c(l2, u1), c("lower", "upper"))) else
      return(setNames(c(l1, u1), c("lower", "upper"))) #by construction: l2<l1, u2<u1
  }
  #
  if (type == "linear") {
    method = match.arg(method)
    if (!(method %in% c("plugin", "boot", "pretest")))
      stop("method for test can only be plugin, boot or pretest")
    if (method == "plugin" | method == "boot") {
      lc.res = lcor.comp(xx, type="linear")
      rho1 = lc.res[1]
      rho2 = lc.res[2]
      lc = lc.res[3]
      n = dim(xx)[1]
      alpha = 1 - conf.level
      if (method == "plugin") Sigma = Sigma.est(xx) #general asymptotic theory
      if (method == "boot") { #bootstrap estimate of covariance matrix
        Sigma = n * cov( boot::boot(xx, statistic = function(xx, id) lcor.comp(xx[id, ], type="linear"), R = R)$t[,1:2] )
      }
      #case |rho1| != |rho2|
      if (abs(rho1) >= abs(rho2)) {
        se = sqrt(Sigma[1,1])
      } else {
        se = sqrt(Sigma[2,2])
      }
      z = qnorm(1 - alpha/2)
      l1 = max( lc - z * se / sqrt(n), 0)
      u1 = min( lc + z * se / sqrt(n), 1)
      #case |rho1|=|rho2|>0
      F = function(x, Sigma, p) {
        sigma1 = sqrt(Sigma[1,1])
        sigma2 = sqrt(Sigma[2,2])
        tau = Sigma[1,2] / ( sigma1 * sigma2 )
        alpha1 = (sigma1/sigma2-tau) / sqrt(1-tau^2)
        alpha2 = (sigma2/sigma1-tau) / sqrt(1-tau^2)
        0.5 * ( sn::psn(x,0,sigma1,alpha1) + sn::psn(x,0,sigma2,alpha2) ) - p
      }
      si = 5 * max(Sigma[1,1],Sigma[2,2]) #guess for search interval
      c.u = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = alpha/2)$root
      c.o = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = 1-alpha/2)$root
      l2 = max( lc - c.o / sqrt(n), 0)
      u2 = min( lc - c.u / sqrt(n), 1)
      if (con==TRUE) return(setNames(c(l2, u1), c("lower", "upper"))) else
        return(setNames(c(l1, u1), c("lower", "upper"))) #by construction: l2<l1, u2<u1
    }
    #
    if (method == "pretest") { #asymptotic theory with pretest
      lc.res = lcor.comp(xx, type="linear")
      rho1 = lc.res[1]
      rho2 = lc.res[2]
      lc = lc.res[3]
      n = dim(xx)[1]
      alpha = 1 - conf.level
      Sigma = Sigma.est(xx) #general asymptotic theory
      #perform pretest for equality of |rho1| and |rho2|
      sigma12 = Sigma[1,2]
      #sign change in asymptotic covariance if -rho1=rho2!=0 (case 2)
      if (sign( rho1*rho2 ) < 0) sigma12 = -sigma12
      S = sqrt(n) * (abs(rho1) - abs(rho2)) / sqrt(Sigma[1,1] - 2*sigma12 + Sigma[2,2])
      if (abs(S) < qnorm(1-alpha/2)) { #assume abs(rho1) = abs(rho2) (!=0)
        F = function(x, Sigma, p) {
          sigma1 = sqrt(Sigma[1,1])
          sigma2 = sqrt(Sigma[2,2])
          tau = Sigma[1,2] / ( sigma1 * sigma2 )
          alpha1 = (sigma1/sigma2-tau) / sqrt(1-tau^2)
          alpha2 = (sigma2/sigma1-tau) / sqrt(1-tau^2)
          0.5 * ( sn::psn(x,0,sigma1,alpha1) + sn::psn(x,0,sigma2,alpha2) ) - p
        }
        si = 5 * max(Sigma[1,1],Sigma[2,2]) #guess for search interval
        c.u = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = alpha/2)$root
        c.o = uniroot(F, c(-si,si), extendInt="upX", Sigma=Sigma, p = 1-alpha/2)$root
        l = max( lc - c.o / sqrt(n), 0)
        u = min( lc - c.u / sqrt(n), 1)
        return( setNames(c(l, u), c("lower", "upper")) )
      }
      #if tests decides for abs(rho1) != abs(rho2)
      if (abs(rho1) >= abs(rho2)) {
        rho = abs(rho1)
        se = sqrt(Sigma[1,1])
      } else {
        rho = abs(rho2)
        se = sqrt(Sigma[2,2])
      }
      z = qnorm(1 - alpha/2)
      l = max( rho - z * se / sqrt(n), 0)
      u = min( rho + z * se / sqrt(n), 1)
      return( setNames(c(l, u), c("lower", "upper")) )
    }
  }
}
