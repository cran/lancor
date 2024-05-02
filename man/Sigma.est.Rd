\name{Sigma.est}
\alias{Sigma.est}
\title{Covariance matrix of components of Lancaster correlation coefficient.}
\description{
Estimate of covariance matrix of the two components of Lancaster correlation.
Lancaster correlation is a bivariate measures of dependence.
}
\usage{
Sigma.est(xx)
}
\arguments{
  \item{xx}{a matrix or data frame with two columns.
}
}
\value{
\verb{Sigma.est} returns the estimated covariance matrix.
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
}
\author{Hajo Holzmann, Bernhard Klar}
\seealso{\code{\link{lcor.ci}}}
\examples{
Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
R <- chol(Sigma)
n <- 1000
x <- matrix(rnorm(n*2), n) %*% R
nu <- 8
y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
Sigma.est(y)
}
