lcor = function(x, y = NULL, type = c("rank", "linear")) {
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
    y = qnorm( (rank(y)-0.5) / length(x) )
    lc = max( abs(cor(x,y)), abs(cor(x^2,y^2)) )
    return(lc)
  }
  if (type == "linear") { 
    x = scale(x)
    y = scale(y)
    lc = max( abs(cor(x,y)), abs( cor(x^2,y^2) ) )
    return(lc)
  }
}
