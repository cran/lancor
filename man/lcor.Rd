\name{lcor}
\alias{lcor}
\title{Lancaster correlation}
\description{
Computes the Lancaster correlation coefficient.
}
\usage{
lcor(x, y = NULL, type = c("rank", "linear"))
}
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
}
  \item{y}{NULL (default) or a vector with same length as x.
}
  \item{type}{a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
}
}
\value{
\verb{lcor} returns the sample Lancaster correlation.
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
}
\author{Hajo Holzmann, Bernhard Klar}
\seealso{\code{\link{lcor.comp}}, \code{\link{lcor.ci}}, \code{\link{lcor.test}}}
\examples{
Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
R <- chol(Sigma)
n <- 1000
x <- matrix(rnorm(n*2), n) %*% R #multivariate normal
lcor(x, type = "rank")
lcor(x, type = "linear")

x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu)
cor(y[,1], y[,2], method = "spearman")
lcor(y, type = "rank")
}
