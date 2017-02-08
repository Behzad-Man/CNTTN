\name{rcpptnm}
\alias{rcpptnm}
\title{Conditional simulation from truncated multivariate normal (TMN) distribution}
\description{
  This function performs conditional univariate simulation from TMN distribution which is required by the slice sampler for
  \eqn{\lambda} and \eqn{U} nodes.
}
\usage{
rcpptnm(sigma, a, mu, sig, ssig, low, high)
}
\arguments{
  \item{sigma}{inverse of the spatial correlation matrix.}
  \item{a}{matrix of latent nodes which have the TMN distribution.}
  \item{mu}{matrix of nodes mean.}
  \item{sig}{vector of diagonal values of the spatial covariance matrix.}
  \item{ssig}{vector of square root of \code{sig}.}
  \item{low}{matrix of lower-truncated values.}
  \item{high}{matrix of upper-truncated values.}
}
\details{
Given the arguments values, this function implements the following R code in C++
	\deqn{for(j in 1:n){for(i in 1:N){a[i,j] <- rtn(mu[i,j] - sig[i]*sigma[i,-i] %*% (a[-i,j] + mu[-i,j]), ssig[j],low[i,j],high[i,j])
	}
	}
	}
and returns \eqn{a}.
}
\value{
  'rcpptnm' returns random deviates with TMN distribution.
}
\author{
  Behzad Mahmoudian \email{b.mahmoudian@qom.ac.ir}
}
