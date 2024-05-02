lcor.comp = function(x, y = NULL, type = c("rank", "linear"), plot = FALSE) {
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
  if (type == "rank") {
    x = qnorm( (rank(x)-0.5) / length(x) )
    y = qnorm( (rank(y)-0.5) / length(y) )
    rho1 = cor(x,y)
    rho2 = cor(x^2,y^2)
    lc = max( abs(rho1), abs(rho2) )
  }
  if (type == "linear") {
    x = scale(x)
    y = scale(y)
    rho1 = cor(x,y)
    rho2 = cor(x^2,y^2)
    lc = max( abs(rho1), abs(rho2) )
  }
  #
  if (plot==TRUE) {
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,2))
    plot(y ~ x, pch=16, xlab="x.standardized", ylab="y.standardized")
    plot(y^2 ~ I(x^2), cex=0.7, pch=16, xlab="square of x.standardized",
         ylab="square of y.standardized")
  }
  return( setNames( c(rho1, rho2, lc), c("rho1", "rho2", "lcor") ) )
}
