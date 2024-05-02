\name{lcor.comp}
\alias{lcor.comp}
\title{Lancaster correlation and its components}
\description{
Computes the Lancaster correlation coefficient and its components.
}
\usage{
lcor.comp(x, y = NULL, type = c("rank", "linear"), plot = FALSE)
}
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
}
  \item{y}{NULL (default) or a vector with same length as x.
}
  \item{type}{a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
}
  \item{plot}{logical; if TRUE, scatterplots of the transformed x and y values and of their squares are drawn.
}
}
\value{
\verb{lcor.comp} returns a vector containing the two components \verb{rho1} and \verb{rho2} and the sample Lancaster correlation.
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
}
\author{Hajo Holzmann, Bernhard Klar}
\seealso{\code{\link{lcor}}, \code{\link{lcor.ci}}, \code{\link{lcor.test}}}
\examples{
Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
R <- chol(Sigma)
n <- 1000
x <- matrix(rnorm(n*2), n) %*% R
nu <- 8
y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
cor(y[,1], y[,2])
lcor.comp(y, type = "linear")

x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
cor(y[,1], y[,2], method = "spearman")
lcor.comp(y, type = "rank", plot = TRUE)
}
