\name{ace.test}
\alias{ace.test}
\title{ACE permutation test of independence}
\description{
Performs a permutation test of independence using ace in package acepack.
ace stands for alternating conditional expectations.
}
\usage{
ace.test(x, y = NULL, nperm = 999)
}
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
}
  \item{y}{NULL (default) or a vector with same length as x.
}
  \item{nperm}{number of permutations.
}
}
\value{
A list containing the following components:
  \item{ace }{the value of the test statistic.}
  \item{pval }{the p-value of the test.}
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
}
\author{Hajo Holzmann, Bernhard Klar}
\seealso{\code{\link{lcor.test}}}
\examples{
n <- 200
x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
cor.test(y[,1], y[,2], method = "spearman")
ace.test(y)
}
