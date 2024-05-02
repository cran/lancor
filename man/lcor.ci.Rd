\name{lcor.ci}
\alias{lcor.ci}
\title{confidence intervals for the Lancaster correlation coefficient
}
\description{
Computes confidence intervals for the Lancaster correlation coefficient.
Lancaster correlation is a bivariate measures of dependence.
}
\usage{
lcor.ci(x, y = NULL, conf.level = 0.95, type = c("rank", "linear"), con = TRUE,
        R = 1000, method = c("plugin", "boot", "pretest"))
}
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
}
  \item{y}{NULL (default) or a vector with same length as x.
}
  \item{conf.level}{confidence level of the interval.
}
  \item{type}{a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
}
  \item{con}{logical; if TRUE (default), conservative asymptotic confidence intervals are computed.
}
  \item{R}{number of bootstrap replications.
}
  \item{method}{a character string indicating how the asymptotic covariance matrix is computed if type ="linear". One of "plugin" (default), "boot" or "symmetric": can be abbreviated.
}
}
\value{
\verb{lcor.ci} returns a vector containing the lower and upper limits of the confidence interval.
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
}
\author{Hajo Holzmann, Bernhard Klar
}
\seealso{\code{\link{lcor}}, \code{\link{lcor.comp}}, \code{\link{lcor.test}}}
\examples{
n <- 1000
x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu) # multivariate t
lcor(y, type = "rank")
lcor.ci(y, type = "rank")
}
