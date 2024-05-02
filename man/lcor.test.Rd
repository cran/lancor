\name{lcor.test}
\alias{lcor.test}
\title{Lancaster correlation test}
\description{
Lancaster correlation test of bivariate independence.
Lancaster correlation is a bivariate measures of dependence.
}
\usage{
lcor.test(x, y = NULL, type = c("rank", "linear"), nperm = 999,
          method = c("permutation", "asymptotic", "symmetric"))
}
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
}
  \item{y}{NULL (default) or a vector with same length as x.
}
  \item{type}{a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
}
  \item{nperm}{number of permutations.
}
  \item{method}{a character string indicating how the p-value is computed if type ="linear". One of "permutation" (default), "asymptotic" or "symmetric": can be abbreviated.
}
}
\value{
  A list containing the following components:
  \item{lcor }{the value of the test statistic.}
  \item{pval }{the p-value of the test.}
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
}
\author{Hajo Holzmann, Bernhard Klar}
\seealso{\code{\link{lcor}}, \code{\link{lcor.comp}}, \code{\link{lcor.ci}}}
\examples{
n <- 200
x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu)
cor.test(y[,1], y[,2], method = "spearman")
lcor.test(y, type = "rank")
}
